#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdint.h>
#include <kseq.h>
#include <hfile.h>
#include <sam.h>
#include <errno.h>
#include "defs.h"
#include "queue.h"
#include "argtable3.h"

KSEQ_INIT(gzFile, gzread);

#define MAX_STEP        5000000
#define MAX_INPUT_SHOW  300
#define MAX_LINE_LENGTH 80

// Externs declared elsewhere.
extern void maxreps(head_t* headnode, char *seq, int length);
extern void mismatchProgram(detected_repeat* dr, char *seq, int length);
extern void printNbOfReps(void);    // printing the number of repetitions found
extern int tooBigReps;              // flag if there exist repetitions which don't fit to the window
extern int actPerBnd;               // this variables are to insure

struct arg_lit *stream, *help, *version;
struct arg_str *level, *file;
struct arg_file *o;
struct arg_end *end;
struct arg_dbl *compression_threshold;
struct arg_int *retains;

/* Parameters - check these for writes! These need to be reset every time the code runs. */
limits lim;
int from;
int to;
int toWasSpecified = NO;
int maxsizeWasSpecified = NO;
int maxperiodWasSpecified = NO;
int start_pstn;                    // window start position
int step = -1;                     // window size
int LastWindow = NO;               // flag indicating whether the window treated is the last one
int maxPer, minPer, dblMinPer, oldPerFactorBound, maxNumErr;
int noprint = 0;                   // flag whether reps themselves should be printed
FILE *output_file;                 // output file (for xml output)
char *output_file_name;
int xmloutput = NO;                // no xml output by default
time_t now;
int allowsmall = NO;    // for indicating small repetitions which may be aleatory

/* Variables */
int inputLength;       // sequence length or file size
int nrRep;            // number of repetitions output
int maxPer, minPer, dblMinPer, oldPerFactorBound;
char nextLetterCheck = '\0', prevLetterCheck = '\0';  // letters immediately following and preceding the current window

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _b : _a; })

// This macro handles all processing.
// Note that maxreps is currently in superstr.c, and mismatchProgram is currently in FndReps.c
#define PROCESS_ACGT(word, length) {                    \
    dblMinPer=2*lim.min_period;                         \
    if ((length)>=dblMinPer) {                            \
        if (lim.err_number == 0)                        \
            maxreps(word,length);                       \
        else                                            \
            mismatchProgram(word,length);               \
        printNbOfReps();                                \
    } else                                              \
        printf("Processed sequence is too short\n\n");  \
}

void init_limits() { // Initialise limits object with sensible defaults for limits (ie, no limits).
    lim.err_number = 0;
    lim.min_period = lim.min_size = 1;
    lim.min_exponent = 1.0f;
    lim.max_period  = lim.max_size = -1;
}

// Wrap up print outputs.
// TODO: Fix this.
void printNbOfReps() {
    // If a repeat was found, print the count line.
    if (nrRep) {
        printf(" ---------------------------------------------------------------------------------------------\n");
    }
    printf("RESULT: Detected %d repeat(s) in the processed sequence.\n\n", nrRep);
    if (tooBigReps) {
        printf("Warning: repeats spanning beyond the window have been detected\n\n");
    }
}

void free_queue(head_t* head) {
    struct s_repeat_queue_node* currnode = NULL;
    while (!(TAILQ_EMPTY(head))) {
        currnode = TAILQ_FIRST(head);
        TAILQ_REMOVE(head, currnode, nodes);
        detected_repeat rpt = currnode->rpt;
        free(rpt.motif);
        free(currnode);
    }
}

/* ************************************************************************* */

// CLI input, start here.
void ProcessSeqFromCommandLine(char *seq, int from, head_t* repeat_head) { // The function handles input.
    // Previous and next letters are null terminators. TODO: work out how to handle better.
    // prevLetterCheck = '\0';
    // nextLetterCheck = '\0';
    start_pstn=from+2; // Start position set by user in from variable
    int final_pstn = (toWasSpecified) ? to : inputLength;
    final_pstn+=2;
    char inputSeq[final_pstn+1];
    // char inputSeq[strlen(seq)+3];
    // Input string size to copy up.
    memcpy(&inputSeq[2], seq, inputLength);
    inputSeq[start_pstn - 1] = inputSeq[start_pstn - 2] = '\0';
    inputSeq[final_pstn] = '\0';
    nrRep = 0;
    tooBigReps = NO;
    LastWindow = YES;
    dblMinPer=2*lim.min_period;
    if (final_pstn-start_pstn>=dblMinPer) {
        if (lim.err_number == 0) {
            maxreps(repeat_head, inputSeq + start_pstn, final_pstn - start_pstn);
        } else {
            detected_repeat app_repeat = {.size=0, .period=0, .errorRate=0.0f};
            mismatchProgram(&app_repeat, inputSeq + start_pstn, final_pstn - start_pstn);
        }
    } else {
        printf("Processed sequence is too short\n\n");
    }
}

void process_string(char* lf_str, int inputLength, head_t* repeat_head) {
    init_limits(); // initialize borderline values
    lim.min_period = 2;
    lim.max_period = 20;
    lim.min_exponent = 3;
    /* 1. Parsing the arguments of the command line */
    allowsmall = NO;
    // For documentation on resolution see superstr documentation at https://mreps.univ-mlv.fr/tutorial.html#resolution
    if (lim.max_period == -1)       // MAXPERIOD was not specified
        lim.max_period = inputLength / 2;
    else {
        maxperiodWasSpecified = YES;
        if (lim.max_period > inputLength / 2) lim.max_period = inputLength / 2;
        if (lim.max_period < lim.min_period) { // TODO: Convert to exception.
            fprintf(stderr, "Error: Maximal period must be greater or equal to the minimal period! \n");
            PROCESS_ERROR(output_file, output_file_name, 4)
            exit(4);
        }
    }
    if (lim.max_size == -1)         // MAXSIZE was not specified
        lim.max_size = inputLength;
    else {
        maxsizeWasSpecified = YES;
        if (lim.max_size < lim.min_size) { // TODO: Convert to exception.
            fprintf(stderr, "Error: Maximal size must be greater than or equal to the minimal size! \n");
            PROCESS_ERROR(output_file, output_file_name, 5)
            exit(5);
        }
    }
    if (from > inputLength) { // TODO: Convert to exception.
        fprintf(stderr, "Error: FROM(=%d) is too big!\n", from);
        PROCESS_ERROR(output_file, output_file_name, 6)
        exit(6);
    }
    if (to > 0) { // Validate TO when specified:
        toWasSpecified = YES;
        if (to < from) { // TODO: Convert to exception.
            fprintf(stderr, "Error: TO(=%d) must be greater or equal than FROM(=%d)\n", to, from);
            PROCESS_ERROR(output_file, output_file_name, 7)
            exit(7);
        }
        if (to > inputLength) { // TODO: Convert to exception.
            printf("Warning: TO(=%d) exceeds sequence length\n", to);
            printf("From: %d\n", from);
            printf("IL: %d\n", inputLength);
            to = inputLength;
        }
    } else { // TO was not specified
        to = inputLength;
    }
    from--; /* Decrement FROM to translate 'natural enumeration' to C indexing */
    if (lim.err_number > 0 && lim.err_number >= inputLength / 2 - 1) { // TODO: Change to exception
        fprintf(stderr, "Error: number of errors (=%d) is too big with respect to the sequence length\n",
                lim.err_number);
        PROCESS_ERROR(output_file, output_file_name, 71)
        exit(71);
    }
    if (lim.min_period < 1) { // Set the min period to 1.
        lim.min_period = 1;
    }
    ProcessSeqFromCommandLine(lf_str, from, repeat_head);
}

int check_ns(char* seq, int seqLen) {
    int n_count = 0;
    char symbol;
    int pos = 0;
    while (pos != seqLen) {
        symbol = seq[pos];
        if (symbol == 'n' || symbol == 'N') {
            n_count++;
            if (1.00 * n_count / (seqLen) > LIMIT_N_PROPORTION) {
                //fprintf(stderr,"\n%dns in %s", n_count, seq);
                return 0;
            }
        }
        pos++;
    }
    return 1;
}
/* The main program */

int process_fastq(char** fastq_files, int filecount, char* outfile, int write_retains, double threshold) {
    // Declare and init variable for counting processed reads.
    long long totalcount = 0;
    // Declare files for sample profile:
    gzFile retain_names;
    gzFile read_str_fp;
    // Build file path for output file:
    char read_str_path[sizeof(outfile)+50];
    strcpy(read_str_path, outfile);
    strcat(read_str_path, "per_read.txt.gz");
    read_str_fp = gzopen(read_str_path, "wb");
    // If write retain flag is on, build the file path for retained read ids, and initialise as required:
    if (write_retains) {
        char file_path[sizeof(outfile) + 50];
        strcpy(file_path, outfile);
        strcat(file_path, "retain_read_ids.txt.gz");
        retain_names = gzopen(file_path, "wb");
    }
    // Declare index for file iteration, then iterate over each file in the input file list.
    int file_idx;
    for (file_idx = 0;file_idx< filecount; file_idx++) {
        gzFile fp = NULL;
        int fno = open(fastq_files[file_idx], O_RDONLY);
        fp = gzdopen(fno, "r");
        // Define pointers to headers/seqs:
        kseq_t *seq;
        seq = kseq_init(fp);
        z_stream defstream;
        defstream.zalloc = Z_NULL;
        defstream.zfree = Z_NULL;
        defstream.opaque = Z_NULL;
        deflateInit(&defstream, Z_BEST_SPEED);
        while (kseq_read(seq) >= 0) {
            defstream.avail_in = (uInt) strlen(seq->seq.s); // size of input, string + terminator
            defstream.next_in = (Bytef *) seq->seq.s; // input char array
            char a[seq->seq.l];
            defstream.avail_out = (uInt) sizeof(a); // size of output
            defstream.next_out = (Bytef *) a; // output char array
            deflate(&defstream, Z_FINISH);
            double uncompressed_size = (double) defstream.total_in;
            double compressed_size = (double) defstream.total_out;
            deflateReset(&defstream);
            double comp_1;
            comp_1 = compressed_size / uncompressed_size;
            totalcount++;
            int n_check_read_1 = check_ns(seq->seq.s, strlen(seq->seq.s));
            if ((n_check_read_1 && (comp_1 < threshold))) {
                if (write_retains) {
                    gzprintf(read_str_fp, "@%s\n%s\n%s\n%s\n", seq->name.s, seq->seq.s, "+", seq->qual.s);
                }
                if ((comp_1 < threshold) && n_check_read_1) {
                    int do_write_flag = 1;
                    head_t repeat_head;
                    TAILQ_INIT(&repeat_head);
                    //0 <= return value < max_int, conversion will be safe unless code run on reads with len > INT_MAX.
                    inputLength = (int) strlen(seq->seq.s);
                    from = 1;
                    to = -1;
                    process_string(seq->seq.s, inputLength, &repeat_head);
                    struct s_repeat_queue_node *e = NULL;
                    TAILQ_FOREACH(e, &repeat_head, nodes) {
                        if (do_write_flag == 1) {
                            gzprintf(read_str_fp, "@%s", seq->name.s);
                            do_write_flag = 0;
                        }
                        gzprintf(read_str_fp, "\t%s:%d:%d:%d:%d:%f", e->rpt.motif, e->rpt.size, e->rpt.period,
                                 e->rpt.start, e->rpt.end, e->rpt.errorRate);
                    }
                    free_queue(&repeat_head);
                    if (do_write_flag == 0) {
                        gzprintf(read_str_fp, "\n");
                    }
                    e = NULL;
                }
            }
        }
        //Clean up kseq and zlib objects.
        deflateEnd(&defstream);
        kseq_destroy(seq);
        //Close the input file.
        gzclose(fp);
        printf("Processing complete of file %d of %d. Closing file handles.\n", file_idx, filecount);
    }
    gzprintf(read_str_fp,"Total %lld\n", totalcount);
    // Close the per-read file.
    gzclose(read_str_fp);
    // Check to see whether the write_retains flag was set, close the file handler if so.
    if (write_retains) {
        gzclose(retain_names);
    }
    printf("Successful run; processed %lld reads.\n", totalcount);
    return 0;
}

int process_stream_fastq(char* outfile, int write_retains, double threshold) {
    long long totalcount = 0;
    char header1[1024];
    char header2[1024];
    char qual1[1024];
    char qual2[1024];
    char comment1[1024];
    char comment2[1024];
    char seq1[1024];
    char seq2[1024];
    FILE *outfp = NULL;
    FILE *outfp2 = NULL;
    if (write_retains) {
        // Declare files for fastq dumping
        char file_path[sizeof(outfile) + 50];
        strcpy(file_path, outfile);
        strcat(file_path, "R1_retain.fa");
        outfp = fopen(file_path, "w");
        strcpy(file_path, outfile);
        strcat(file_path, "R2_retain.fa");
        outfp2 = fopen(file_path, "w");
    }
    // Declare file for per-read STR dump:
    gzFile read_str_fp;
    char read_str_path[sizeof(outfile)+50];
    strcpy(read_str_path, outfile);
    strcat(read_str_path, "per_read.txt.gz");
    read_str_fp = gzopen(read_str_path, "wb");
    //// Define pointers to headers/seqs
    z_stream defstream;
    defstream.zalloc = Z_NULL;
    defstream.zfree = Z_NULL;
    defstream.opaque = Z_NULL;
    deflateInit(&defstream, Z_BEST_SPEED);
    while (4 == fscanf(stdin, "%[^\n]\n%[^\n]\n%[^\n]\n%[^\n]\n", header1, seq1, comment1, qual1)) {
        //// zlib struct
        defstream.avail_in = (uInt) strlen(seq1); // size of input, string + terminator
        defstream.next_in = (Bytef *) seq1; // input char array
        char a[strlen(seq1)];
        defstream.avail_out = (uInt) sizeof(a); // size of output
        defstream.next_out = (Bytef *) a; // output char array
        deflate(&defstream, Z_FINISH);
        double uncompressed_size = (double) defstream.total_in;
        double compressed_size = (double) defstream.total_out;
        deflateReset(&defstream);
        double comp_1;
        comp_1 = compressed_size / uncompressed_size;
        totalcount++;
        int n_check_read_1 = check_ns(seq1, strlen(seq1));
        if ((n_check_read_1 && (comp_1 < threshold))) {
            if (write_retains) {
                fprintf(outfp, "@%s\n%s\n%s\n%s\n", header1, seq1, "+", qual1);
            }
            if ((comp_1 < threshold) && n_check_read_1) {
                int do_write_header_flag = 1;
                head_t repeat_head;
                TAILQ_INIT(&repeat_head);
                //0 <= return value < max_int, conversion will be safe unless code run on reads with len > INT_MAX.
                inputLength = (int) strlen(seq1);
                from = 1;
                to = -1;
                process_string(seq1, inputLength, &repeat_head);
                struct s_repeat_queue_node *e = NULL;
                TAILQ_FOREACH(e, &repeat_head, nodes) {
                    if (do_write_header_flag == 1) {
                        gzprintf(read_str_fp, "@%s", header1);
                        do_write_header_flag = 0;
                    }
                    gzprintf(read_str_fp, "\t%s:%d:%d:%d:%d:%f", e->rpt.motif, e->rpt.size, e->rpt.period,
                             e->rpt.start, e->rpt.end, e->rpt.errorRate);
                }
                free_queue(&repeat_head);
                if (do_write_header_flag == 0) {
                    gzprintf(read_str_fp, "\n");
                }
                e = NULL;
            }
        }
    }
    deflateEnd(&defstream);
    //Check to see whether the write_retains flag was set, close the file handler if so.
    if (write_retains) {
        fclose(outfp);
        fclose(outfp2);
    }
    gzclose(read_str_fp);
    printf("Successful run; processed %lld reads.\n", totalcount);
    return 0;
}

int process_bam_sam(char* bamfile, char* outfile, int write_retains, double threshold) {
    long long totalcount = 0;
    FILE *outfp;
    //Note: htslib's hts_open handles BAM, SAM and CRAM
    samFile *fp_in = hts_open(bamfile,"r"); //open BAM/CRAM file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    if (write_retains) {
        // Declare files for fastq dumping
        char file_path[sizeof(outfile) + 50];
        strcpy(file_path, outfile);
        strcat(file_path, "retain_reads.txt");
        outfp = fopen(file_path, "w");
    }
    // Declare file for per-read STR dump:
    gzFile read_str_fp;
    char read_str_path[sizeof(outfile)+50];
    strcpy(read_str_path, outfile);
    strcat(read_str_path, "per_read.txt.gz");
    read_str_fp = gzopen(read_str_path, "wb");
    z_stream defstream;
    defstream.zalloc = Z_NULL;
    defstream.zfree = Z_NULL;
    defstream.opaque = Z_NULL;
    deflateInit(&defstream, Z_BEST_SPEED);
    int r;
    while((r = sam_read1(fp_in,bamHdr,aln)) >= 0){
        // int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        // char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        char *name = bam_get_qname(aln);
        // Get read length
        uint32_t len = aln->core.l_qseq; //length of the read.
        uint8_t *q = bam_get_seq(aln); //encoded sequence
        size_t flag = aln->core.flag;
        // Must be not an optical/PCR duplicate, a secondary or supplementary alignment. This should leave unmapped and
        // primary alignments only.
        if ((flag & BAM_FSECONDARY) || (flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }
        char *qseq = (char *)malloc(len+1);
        for(int i=0; i < len; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        defstream.avail_in = len; // size of input, string + terminator
        defstream.next_in = (Bytef *) qseq; // input char array
        char a[len];
        defstream.avail_out = (uInt) sizeof(a); // size of output
        defstream.next_out = (Bytef *) a; // output char array
        deflate(&defstream, Z_FINISH);
        double uncompressed_size = (double) defstream.total_in;
        double compressed_size = (double) defstream.total_out;
        deflateReset(&defstream);
        double comp_1;
        comp_1 = compressed_size / uncompressed_size;
        totalcount++;
        if (comp_1 < threshold) {
            int n_check_read_1 = check_ns(qseq, len);
            if ((comp_1 < threshold) && n_check_read_1) {
                if (write_retains) {
                    fprintf(outfp, "@%s\n%s\n", name, qseq);
                }
                head_t repeat_head;
                TAILQ_INIT(&repeat_head);
                //0 <= return value < max_int, conversion will be safe unless code run on reads with len > INT_MAX.
                inputLength = (int) len;
                from = 1;
                to = -1;
                process_string(qseq, inputLength, &repeat_head);
                struct s_repeat_queue_node *e = NULL;
                int do_write_header_flag = 1;
                TAILQ_FOREACH(e, &repeat_head, nodes) {
                    if (do_write_header_flag == 1) {
                        gzprintf(read_str_fp, "@%s", name);
                        do_write_header_flag = 0;
                    }
                    gzprintf(read_str_fp, "\t%s:%d:%d:%d:%d:%f", e->rpt.motif, e->rpt.size, e->rpt.period,
                             e->rpt.start, e->rpt.end, e->rpt.errorRate);
                }
                free_queue(&repeat_head);
                if (do_write_header_flag == 0) {
                    gzprintf(read_str_fp, "\n");
                }
                e = NULL;
            }
        }
        free(qseq);
    }
    gzprintf(read_str_fp,"Total %lld\n", totalcount);
    gzclose(read_str_fp);
    bam_destroy1(aln);
    sam_close(fp_in);
    printf("Successful run; processed %lld reads.\n", totalcount);
    return 0;
}

int main(int argc, char *argv[]) {
    void *argtable[] = {
            help = arg_litn(NULL, "help", 0, 1, "display this help and exit"),
            version = arg_litn(NULL, "version", 0, 1, "display version info and exit"),
            level = arg_strn(NULL, "mode", "<fastq|bam>", 0, 1, "type of input data"),
            stream = arg_litn("s", "stream", 0, 1, "run on named streams, not files (see manual for instructions)"),
            o = arg_filen("o", NULL, "myfile", 0, 1, "output directory"),
            file = arg_strn(NULL, NULL, "<file>", 0, 100, "input files (or names of pipes in stream mode)"),
            compression_threshold = arg_dbln("t","threshold","float",0,1,"compression threshold for processing strings"),
            retains = arg_int0("r", "retain", NULL, "write retained reads"),
            end = arg_end(20),
    };
    int exitcode = 0;
    char progname[] = "superstr";
    int nerrors;
    level->sval[0] = "";
    compression_threshold->dval[0] = 0.40;
    retains->ival[0] = 0;
    nerrors = arg_parse(argc,argv,argtable);
    // Show help
    if (help->count > 0) {
        printf("Usage: %s\n", progname);
        arg_print_syntax(stdout, argtable, "\n");
        printf("Rapid STR characterisation in NGS data.\n");
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        exitcode = 0;
        goto exit;
    }
    // Show error text for arguments - update to be a bit more sophisticated.
    if (nerrors > 0) {
        arg_print_errors(stdout, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        exitcode = 1;
        goto exit;
    }
    char** file_set = (char **) file->sval;
    char* output_path = (char *) o->filename[0];
    char* mode = (char *) level->sval[0];;
    for (int idx = 0; idx < sizeof(file_set); idx++) {
        if (access( file_set[idx], F_OK ) == -1 ) {
            if (file_set[idx][0] != '\0')  {
                printf("ERROR: Input file %s does not appear to exist.\n", file_set[idx]);
                exitcode = 1;
                goto exit;
            }
        }
    }
    if (access( output_path, F_OK ) == -1 ) {
        printf("ERROR: Output file %s does not appear to exist.\n", output_path);
        exitcode = 1;
        goto exit;
    }
    if (stream->count == 0) {
        if ((strcmp(mode, "bam") == 0) && (strcmp(mode, "fastq") == 0)) {
            printf("ERROR: Mode must be one of 'bam' or 'fastq', got '%s'\n", mode);
            exitcode = 1;
            goto exit;
        }
    }
    if (stream->count > 0) {
        exitcode = process_stream_fastq(output_path, retains->ival[0], compression_threshold->dval[0]);
    } else if (strcmp(mode, "bam") == 0) {
        exitcode = process_bam_sam(file_set[0], output_path, retains->ival[0], compression_threshold->dval[0]);
    } else {
        exitcode = process_fastq(file_set, file->count, output_path, retains->ival[0], compression_threshold->dval[0]);
    }
    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return exitcode;
}


