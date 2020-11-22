# Processing samples with superSTR

superSTR follows a simple process for processing sequencing data.

For each read within your FASTQ/BAM file, superSTR:
1. (BAM only) Excludes any read with [SAM flag](https://samtools.github.io/hts-specs/) 256 ('secondary alignment') or 2048 ('supplementary alignment') set.
2. Compresses the read sequence with zlib.
3. Calculates a ratio of the in-memory size of the compressed read sequence to the size of the uncompressed read sequence.
4. If that ratio is less than the user-specified threshold, it runs mreps on the read with the default mreps parameters specified below.

## Threshold selection

The main difficulty in running superSTR is in identifying the appropriate compression ratio to use as a threshold. This ratio is dependent on the read length and the GC content of the sequence of interest, and must be 

The recommended thresholds for maximising accuracy and precision/recall for human genomic data (42% GC content) at common read lengths are:

| Read Length | Accuracy | Precision + Recall |
|-------------|----------|--------------------|
| 75 | 0.64 | 0.64 |
| 150 | 0.50 | 0.49 |
| 200 | 0.46 | 0.46 |


## Running superSTR

**FASTQ:** `superstr --mode=fastq -o output_dir/ -t 0.64 input_1.fastq.gz input_2.fastq.gz`

**BAM:** `superstr --mode=bam -o output_dir/ -t 0.64 input.bam`

Note: The -o flag is an output prefix rather than a file path; "per_read.txt.gz" is appended to the prefix. 

This can produce slightly counterintuitive behaviour:

`-o output` will produce a file called `outputper_read.txt.gz`

`-o output/` will produce a file called `per_read.txt.gz` in the `output/` directory

## Advanced or non-standard usage:

The [threshold documentation](docs/THRESHOLD.md) contains a full guide to superSTR thresholds, including an explanation of this table, how it was generated, and information for non-human and general thresholds.

The advanced execution documentation will contain a full list of superSTR input and output modes, including instructions for running superSTR on pipes and named pipes (for example with `fastq-dump` from the SRA toolkit or the outputof the `gsutil cat` command on Google Cloud).
