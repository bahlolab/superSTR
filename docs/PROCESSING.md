# Processing samples with superSTR

superSTR follows a simple process for processing sequencing data.

For each read within your FASTQ/BAM/CRAM file, superSTR:
1. (BAM/CRAM only) Excludes any read with [SAM flag](https://samtools.github.io/hts-specs/) 256 ('secondary alignment') or 2048 ('supplementary alignment') set.
2. Compresses the read sequence with zlib.
3. Calculates a ratio of the in-memory size of the compressed read sequence to the size of the uncompressed read sequence.
4. If that ratio is less than the user-specified threshold, it runs mreps on the read with the default mreps parameters specified below.

## Threshold selection

The main difficulty in running superSTR is in identifying the appropriate compression ratio to use as a threshold. This ratio is dependent on the read length and the GC content of the sequence of interest, and must be 

Thresholds for maximising accuracy in some common genetic contexts at common read lengths are:

| Read Length | H. sapiens (~41% GC) | S. coelicolor (72% GC) | P. falciparum (13.5% GC) | 50% GC | 
|-------------|----------|----------|----------|----------|
| 75 | 0.63 | 0.63 | 0.61 | 0.63 |
| 100 | 0.56 | 0.55 | 0.54 | 0.56 |
| 125 | 0.52 | 0.51 | 0.49 | 0.52 |
| 150 | 0.49 | 0.48 | 0.47 | 0.49 |
| 200 | 0.46 | 0.45 | 0.43 | 0.46 |
| 250 | 0.44 | 0.42 | 0.40 | 0.44 |
| 300 | 0.42 | 0.41 | 0.39 | 0.42 |

You can generate similar tables for arbitrary read lengths (and for maximising precision/recall) using simulation.py in the Python script folder.

Tuning these thresholds down may significantly improve run times at the expense of short repeat sequences. An extension of the simulation method to assist in estimating the impact of tuning these thresholds is in development.

## Running superSTR

**FASTQ:** `superstr --mode=fastq -o output_dir/ -t 0.64 input_1.fastq.gz input_2.fastq.gz`

**BAM:** `superstr --mode=bam -o output_dir/ -t 0.64 input.bam`

**CRAM:** `superstr --mode=bam -o output_dir/ -t 0.64 input.cram`

superSTR should respect environment settings of REF_PATH and REF_CACHE as outlined in the [htslib documentation]([http://www.htslib.org/workflow/](https://www.htslib.org/workflow/cram.html)). Execution times can vary for CRAM if these settings aren't used due to the need to download and cache sequences from the EBI servers.

Note: The -o flag is an output prefix rather than a file path; "per_read.txt.gz" is appended to the prefix. The current version of superSTR does not create output directories if they are not found.

This can produce slightly counterintuitive behaviour:

`-o output` will produce a file called `outputper_read.txt.gz`

`-o output/` will produce a file called `per_read.txt.gz` in the `output/` directory


## Advanced or non-standard usage:

superSTR can also be run on named pipes (or [FIFO](https://man7.org/linux/man-pages/man7/fifo.7.html)s), for example with `fastq-dump` from the SRA toolkit or the output of the `gsutil cat` command on Google Cloud. If using this functionality it's recommended that you treat paired fastq files as single-end fastq and merge them at the end of analysis.
