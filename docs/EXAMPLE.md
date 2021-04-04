**This file is being worked on; there may be errors in this example. This warning will be removed once the example is completed and tested fully.**

# Example workflow - SCA3 RNA-seq analysis

This example workflow steps through all core functions of processing a public SCA3 RNA-seq dataset using superSTR from the NCBI Sequence read archive, as outlined in the paper. The data is 24 bulk RNA-seq experiments performed on peripheral blood monocyctes from 12 individuals with SCA3 and 12 matched controls from a study of spinocerebellar ataxia type 3 in patients from mainland China (from [Li, T. et al. Front. Genet., 13 June 2019](https://doi.org/10.3389/fgene.2019.00566)).

# Introduction and requirements

The workflow is written to automate as much of this process as possible on a server running a standard Linux operating system. It has been tested under Centos 8. Running this workload will download a large amount of data (>103GB), then convert it to compressed text.

You'll need about 300 GB of storage available to be able to complete this successfully:

* ~120GB for SRA downloads (.sra and reference files for sra-toolkit)
*

You'll need to make sure the following commands are on your PATH variable. cat, xargs and gzip should be already installed; you can check each command by running the command in square brackets after each bullet point):

* prefetch (from sra-toolkit) [prefetch --help]
* fasterq-dump (from sra-toolkit) [fastq-dump --help]
* cat [cat --help]
* xargs [xargs --help]
* gzip [gzip --help]

The below bash scripts are also provided as a Nextflow pipeline to run locally under the same assumptions; you can extend this to use executors (e.g. Slurm, PBS/Torque) as appropriate.

# Overview - superSTR processing

You'll need to download and install superSTR as per the instructions [here](docs/INSTALL.md).

We will walk through the following steps (which are generally applicable to most analyses):

1. [Download data from the SRA](#downloading-data-from-the-sequence-read-archive)
2. [Process reads with superSTR](#process-reads-with-superSTR)
3. Build a manifest file
4. Summarise the superSTR run
5. Perform dataset-level motif screening
6. Perform sample-level outlier detection (optional)
7. Visualise motif screening outputs

# Downloading data from the Sequence Read Archive

**Entry point:** Start in the `examples/` directory of the superSTR git repository.

This process might take some time, depending on your internet connection speed - it downloads ~120GB of files from the SRA, and converts them to gzipped FASTQ files. Make sure that you delete your SRA cache when you're done!

1. Download data:
`prefetch --option-file files/SRP168964_accessions.txt
`

2. Convert data to fastq.gz files:

  Use one of:

  **fastq-dump** - slower, more space efficient:

  `cat files/SRP168964_accessions.txt | xargs -n 1 -P 12 -I % fastq-dump % --split-3 --skip-technical --gzip`

  **fasterq-dump** - this may be faster depending on your configuration, but uses additional space:

  `cat files/SRP168964_accessions.txt | xargs -n 1 -P 1 fasterq-dump`

If you are working with data that needed trimming to uniform read length (or other preprocessing of fastq data), do this here; the example data doesn't need this.

# Process reads with superSTR

**Input:** 48 fastq.gz files in the `examples/` directory, two files per SRR run ID.

**Output:** superSTR processed reads in the `superSTR_rp/` directory; one subdirectory per sample containing a per_read.txt.gz file.

1. Make a directory structure for superSTR to write results to - a root directory, then one subdirectory per sample:

  `mkdir superSTR_rp/`

  `cat files/SRP168964_accessions.txt | xargs -n 1 -I % mkdir superSTR_rp/%/`

2. Run superSTR read processing (the threshold is set using `-t 0.56` for H. sapiens at 100nt read length):

  `cat files/SRP168964_accessions.txt | xargs -n 1 -P 12 -I % superSTR/C/superstr --mode=fastq -o superSTR_rp/%/ -t 0.56 %_1.fastq.gz %_2.fastq.gz`

# superSTR summarisation

**Input:** The output of the previous step.

**Output:** Summarised superSTR data in the `superSTR_processed/` directory

1. Create a manifest file from this experiment's template with correct paths on your machine (you can also edit this manually):

  `pwd | xargs -I % sed 's:xxxxxx:%:g' files/SRP168964_manifest_template.tsv > SRP168964_manifest.tsv`

2. Make a directory for the output of the processing step:

  `mkdir superSTR_processed/`

3. Process the outputs. You can change the number of parallel processes by modifying the `-@ 12` argument:

  `python ../Python/multiparse.py -@ ${N_CORES} --input superSTR_rp/ --output superSTR_processed/ -r`
# superSTR motif screening

# Motif screening

# Outlier detection

# Visualisation
