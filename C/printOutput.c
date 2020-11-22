#include "defs.h"
#include "stdio.h"
#include "string.h"
#include "queue.h"

#define EXPERIMENT
extern int start_pstn, nrRep;
extern char *seq_original;
extern int xmloutput;
extern FILE *output_file;


// DNA Base table from Alex Reynold: https://gist.github.com/alexpreynolds/4f75cab4350e9d937f4a#file-rc-c
static const unsigned char basemap[256] = {
        0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
        16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
        32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
        48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

char *str_reverse_in_place(char *str, int len)
{
    char *p1 = str;
    char *p2 = str + len - 1;
    while (p1 < p2) {
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
    return str;
}

char *str_complement(char *str, char *dest, int len) {
    for (int idx = 0; idx < len; idx++) {
        dest[idx] = basemap[(int)str[idx]];
    }
    return dest;
}
//
//#ifdef EXPERIMENT
//void print_score(int rinitpos, int rendpos, int rlength, int rperiod, float rscore)
///* printing an exact or approximate rep */
//{
//    printf("\nRinitpos: %d, Rendpos: %d, Rperiod: %d\n", rinitpos, rendpos, rperiod);
//    if ((nrRep++)==0 && xmloutput==NO)
//    {
//        printf("   from   ->       to  : \t size\t <per.>\t [exp.]\t\t err-rate \tsequence\n") ;
//        printf(" ---------------------------------------------------------------------------------------------\n");
//    }
//
//    if (xmloutput==YES)
//        fprintf (output_file,
//                 "\t<repeat>\n\t\t<start>%d</start>\n\t\t<end>%d</end>\n\t\t<length>%d</length>\n\t\t<period>%d</period>\n\t\t<exponent>%.2f</exponent>\n\t\t<score>%.3f</score>\n",
//                 start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod,((float) rlength) / rperiod, rscore) ;
//    else
//        printf("%8d  ->  %8d :   \t %d \t <%d> \t [%.2f] \t %.3f \t\t",
//               start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod, ((float) rlength) / rperiod, rscore
//        ) ;
//
//        int per_start,pos_in_per;
///*       int re = (rendpos<rinitpos + MAXDISPLAY)? rendpos : rinitpos+(rperiod-1); */
//        if (xmloutput==YES)
//            fprintf(output_file,"\t\t<sequence>\n");
//
//        for (per_start = rinitpos; per_start <= rendpos; per_start+=rperiod)
//        {
//            if (xmloutput==YES)
//                fprintf(output_file,"\t\t\t<unit>");
//            for (pos_in_per=0; pos_in_per<rperiod && per_start+pos_in_per<=rendpos; pos_in_per++)
//            {
//                if (xmloutput==YES)
//                    fputc(seq_original[per_start+pos_in_per-1],output_file);
//                else
//                    printf ("%c",seq_original[per_start+pos_in_per-1]);
//            }
//            if (xmloutput==YES)
//                fprintf(output_file,"</unit>\n");
//            else
//                printf(" ");
//        }
//        if (xmloutput==YES)
//            fprintf(output_file,"\t\t</sequence>\n");
//
//    if (xmloutput==YES)
//        fprintf(output_file,"\t</repeat>\n");
//    else
//        printf ("\n");
//}
//#endif


#ifdef EXPERIMENT
void print_score(head_t* repeat_head, int rinitpos, int rendpos, int rlength, int rperiod, float rscore) { /* printing an exact or approximate rep */
    char minimal[rperiod+1];
    bzero(minimal, rperiod+1);
    int per_start, pos_in_per;
    // Print repeat split up into repeating motif:
    for (per_start = rinitpos; per_start <= rinitpos+rperiod; per_start += 1) {
        char test_string[rperiod+1];
        bzero(test_string, rperiod+1);
        strncpy(test_string, seq_original + (per_start - 1), rperiod);
        test_string[rperiod] = '\0';
        if (minimal[0] == '\0')
            strcpy(minimal, test_string);
        if (strcmp(minimal, test_string) > 0)
            strcpy(minimal, test_string);
        char complement[rperiod + 1];
        complement[rperiod] = '\0';
        str_complement(test_string, complement, rperiod);
        // TO-DO: Add a flag for paranoia where it's assumed that nucleotides doesn't run 5-3' in the presented string
        //if (strcmp(minimal, complement) > 0) {
        //    strcpy(minimal, complement);
        //}
        //str_reverse_in_place(test_string, rperiod);
        //if (strcmp(minimal, test_string) > 0)
        //    strcpy(minimal, test_string);
        str_reverse_in_place(complement, rperiod);
        if (strcmp(minimal, complement) > 0) {
            strcpy(minimal, complement);
        }
    }
    // At this point we have:
    //  minimal = lex order minimal motif
    //  rperiod = period
    //  rlength = length
    //  rscore = error score
    // Put this into another struct to return.
    struct s_repeat_queue_node* insert_node;
    insert_node = malloc(sizeof(struct s_repeat_queue_node));
    detected_repeat dr;
    dr.period=rperiod;
    dr.motif=strdup(minimal);
    //insert_node->rpt.motif=&minimal[0];
    dr.size=rlength;
    dr.errorRate=rscore;
    dr.start=rinitpos;
    dr.end=rendpos;
    insert_node->rpt=dr;
    TAILQ_INSERT_TAIL(repeat_head, insert_node, nodes);
}
#endif
