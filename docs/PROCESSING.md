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


## 
**FASTQ**


## Advanced or non-standard usage:

The [threshold documentation](docs/THRESHOLD.md) contains a full guide to superSTR thresholds, including an explanation of this table, how it was generated, and information for non-human and general thresholds.
