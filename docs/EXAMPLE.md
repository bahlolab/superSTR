# Example workflow - SCA3 RNA-seq analysis

This example workflow steps through all core functions of processing a sample of a public SCA3 RNA-seq dataset from raw sequence (in gzipped FASTQ format) using superSTR from the NCBI Sequence read archive.

The data is drawn from bulk RNA-seq experiments performed on peripheral blood monocyctes of individuals with SCA3. We take three individuals with SCA3 and three controls from a study of SCA3 in patients from mainland China (from [Li, T. et al. Front. Genet., 13 June 2019](https://doi.org/10.3389/fgene.2019.00566)).

# Introduction and requirements

The workflow is written to automate as much of this process as possible on a server running a standard Linux operating system. It has been tested under Centos 8.

We use a small subset of data from a larger study to reduce the download and compute requirements involved in this demo. Downloading data and converting to FASTQ represents the bulk of the time involved in this introduction.

You'll need about 70 GB of storage available to be able to complete this successfully:

* ~28GB for SRA downloads (.sra and reference files used by the sra-toolkit)
* ~37GB for the tutorial (37GB of fastq.gz files; <5MB for all other files)

You'll need to make sure the following commands are on your PATH variable. cat, xargs and gzip should be already installed; you can check each command by running the command in square brackets after each bullet point):

* prefetch (from sra-toolkit) [`prefetch --help`]
* fastq-dump (from sra-toolkit) [`fastq-dump --help`]
* cat [`cat --help`]
* xargs [`xargs --help`]
* gzip [`gzip --help`]

The below bash scripts will also be provided as a Nextflow pipeline to run locally under the same assumptions soon; you will be able to extend this to use executors (e.g. Slurm, PBS/Torque) as appropriate.

# Overview - superSTR processing

You'll need to download and install superSTR as per the instructions [here](INSTALL.md).

We will walk through the following steps (which are generally applicable to most analyses):

1. [Download data from the SRA](#downloading-data-from-the-sequence-read-archive)
2. [Process reads with superSTR](#process-reads-with-superSTR)
3. [Build a manifest file](#build-a-manifest-file)
4. [Summarise the superSTR run](#summarise-the-run)
5. [Perform motif screening](#perform-motif-screening)
6. [Perform sample-level outlier detection (optional)](#perform-outlier-detection)
7. [Visualise motif screening outputs](#visualise-outputs)

All steps in this tutorial example require that you start in a command prompt with the `examples/` directory of the cloned superSTR git repository as your working directory.

# Downloading data from the Sequence Read Archive

**Generates:** FASTQ files for analysis.

**Notes:** The longest step in this example workflow is that of obtaining data to analyse. We have selected three cases and three controls from a larger experiment to demonstrate superSTR without requiring large amounts of storage space. Samples were chosen as the first three cases and first three controls in the lexicographic order of accessions.

Downloading ~37GB data from the SRA and converting it to fastq files may still take some time depending on your internet connection speed. Make sure that you delete your SRA cache when you're done!

If you were working with data that needed trimming to uniform read length (or other preprocessing of fastq data), do this here; the example data doesn't need this.

**Commands:**

1. Download data:

  `prefetch --option-file files/SRP168964_accessions.txt`

2. Convert data to fastq.gz files.

  *Either*: `fastq-dump` - slower, more space efficient:

  `cat files/SRP168964_accessions.txt | xargs -n 1 -I % fastq-dump % --split-3 --skip-technical --gzip`

  *Or*: `fasterq-dump` - this may be faster depending on your configuration, but uses additional space:

  `cat files/SRP168964_accessions.txt | xargs -n 1 fasterq-dump`


# Process reads with superSTR

**Requires:** FASTQ files for analysis.

**Generates:** superSTR processed reads in the `superSTR_rp/` directory; one subdirectory per sample containing a per_read.txt.gz file.

**Notes:** This is the normal entry point for superSTR analysis. Each sample will take somewhere around 25-35 minutes to run. It is strongly recommended that you run these in parallel on a HPC system; this reduces the total runtime for this step to <40 minutes.

Running superSTR in parallel on a single server may result in bottlenecks in reading data from disk or throttling by user CPU limits on shared machines. You can trial parallel execution on your system by increasing `-P 1` to `-P 6` (or higher). A good way to check for bottlenecks is to monitor the output of the `top` command. superSTR processes *should* be reported as running at 100% CPU.

**Commands:**

1. Make a directory structure for superSTR to write results to - a root directory, then one subdirectory per sample:

  `mkdir superSTR_rp/`

  `cat files/SRP168964_accessions.txt | xargs -n 1 -I % mkdir superSTR_rp/%/`

2. *Either*: Run superSTR read processing one at a time for each sample. The threshold is set using `-t 0.56` for H. sapiens at 100nt read length:

  `cat files/SRP168964_accessions.txt | xargs -n 1 -P 1 -I % superSTR/C/superstr --mode=fastq -o superSTR_rp/%/ -t 0.56 %_1.fastq.gz %_2.fastq.gz`

  *Or*: If you are running on a HPC with Slurm, you may want to run superSTR in parallel using a command such as:

  `cat files/SRP168964_accessions.txt | xargs -n 1 -P 1 -I % sbatch -n1 --mem=2Gb --wrap="superSTR/C/superstr --mode=fastq -o superSTR_rp/%/ -t 0.56 %_1.fastq.gz %_2.fastq.gz"`


# Build a manifest file

**Note:** Manifest files are detailed [here](https://github.com/bahlolab/superSTR/blob/main/docs/POSTPROC.md). This step takes a pre-filled example ([`examples/files/SRP168964_manifest_template.txt`](../examples/files/SRP168964_manifest_template.txt)) for the SCA3 RNA-seq data and sets up the paths to the outputs you've produced.

  You could also do this manually by replacing "xxxxxx" in that file with the relevant path from your machine.

**Commands:**

  `pwd | xargs -I % sed 's:xxxxxx:%/superSTR_rp:g' files/SRP168964_manifest_template.txt > SRP168964_manifest.tsv`

# Summarise the run

**Requires:** raw superSTR data in the `superSTR_rp` directory and a manifest file from the preceding steps.

**Generates:** summarised superSTR data in the `superSTR_processed/` directory.

**Note:** If there is an error, you may need to add the `--clobber` flag to `multiparse.py`, but be careful - this will overwrite existing files in the `--output` path.

**Commands:**

1. Make a directory for the output of the processing step:

  `mkdir superSTR_processed/`

2. Process the outputs. You can change the number of parallel processes by modifying the `-@ 12` argument:

  `python superSTR/Python/multiparse.py -@ 12 --input superSTR_rp/ --output superSTR_processed/ -r 100 -m SRP168964_manifest.tsv`

# Perform motif screening

**Requires:** Summarised superSTR output in the `superSTR_processed/` directory as described in [the preceding step](#summarise-the-run).

**Generates:** A TSV file containing the results of motif screening in `screening.tsv`.

**Notes:** Usually we would only test groups with more than a specified number of samples using the `--min_thresh` argument (this defaults to 3). We have to override this for this tutorial example.

The permutation p-value *may* show a zero value in outputs for this example. This is because the screening code uses [mpmath defaults](https://mpmath.org/doc/current/technical.html), which is 53 bit precision (or decimal precision of 15.95 digits). If the p-value is < 10<sup>-16</sup> then it may display as zero.

The motif screener has help documentation viewable with `python superSTR/Python/screen.py --help`.

**Commands:**

1. Run the motif screener script:

  `python superSTR/Python/screen.py -m SRP168964_manifest.tsv -i superSTR_processed/ -o screening.tsv --controllab Control --swaplab --min_thresh 1`

# Perform outlier detection

**Requires:** Summarised superSTR output in the `superSTR_processed/` directory as described in [the preceding step](#summarise-the-run).

**Generates:** A TSV file containing the results of motif screening in `outliers.tsv`.

**Notes:** The outlier detector has help documentation viewable with `python superSTR/Python/outliers.py --help`

**Commands:**

1. Run the outlier detector script:

  `python superSTR/Python/outliers.py -i superSTR_processed/ -o outliers.tsv --max_motif 6 --bootstrapCI --controllab Control -m SRP168964_manifest.tsv --min_len 75 --max_len 100 -is`

# Visualise outputs

**Requires:** Summarised superSTR output in the `superSTR_processed/` directory as described in [the preceding step](#summarise-the-run).

**Generates:** Split violin plots as shown in Figure 3 of the superSTR manuscript.

**Notes:** Many arguments are being passed directly to matplotlib; they should work as though they were running in matplotlib.

You can specify filetypes for outputs using the extension of the `-o` argument; supported types are [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html). Height and width are in inches, and are passed directly to [figsize](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html).

The visualisation code has help documentation viewable with `python superSTR/Python/Visualise.py --help`

**Commands:**

1. Run the visualisation script on the motif of interest:

  `python superSTR/Python/visualise.py -i superSTR_processed/motifs/3mers/AGC.csv -o AGC.pdf -m SRP168964_manifest.tsv --aff_lab SCA3 --ctrl_lab Control`
