# Installation:

## Requirements

The following instructions have been tested on:

* Debian 10.8 "buster"
* Ubuntu 20.04 LTS "Focal Fossa"
* Centos 8
* Apple OSX 10.14.6 "Mojave" 

superSTR does not currently support Microsoft Windows (there's extensive localisation work needed). If you are running Windows, [a containerised version](#installation-via-docker) of superSTR is available but is untested on that platform. 

## Conda Installation

The easiest way is to install the superSTR conda package

```
conda install -c bioconda -c conda-forge superstr
```

or if you want to create a new environment

```
conda create -n myenv -c bioconda -c conda-forge superstr
conda activate myenv
```

This will install superSTR and all necessary dependencies.

The superSTR binary and the four provided python scripts are then globally available:

```bash
superstr
superstr-multiparse.py
superstr-outliers.py
superstr-screen.py
superstr-visualise.py
```


## Manual installation

Please make sure the following are installed and available on your system prior to installing superSTR:

* [htslib >1.9](https://github.com/samtools/htslib)
* [zlib](https://zlib.net/)
* [Cmake](https://cmake.org/install/) (>3.10)
* [curl](https://curl.se/download.html) and [libcurl](https://curl.se/libcurl/)
* [gcc](https://gcc.gnu.org/install/)
* [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

### read processing

An example installation script is available [in the examples directory](../examples/install_script.sh). The script creates a superSTR directory containing an htslib 1.11 installation, along with a copy of superSTR, then links and compiles the software.

Alternatively:
1) Clone this repo to your machine using `git clone https://github.com/bahlolab/superSTR`.
2) Change into the superSTR directory with `cd superSTR`
3) Set the the HTSLIB_ROOT environment variable with `export HTSLIB_ROOT=<path_to_your_htslib_installation>`
4) Change into the C source directory with `cd C`
5) Run `cmake . && make`
6) Check that superSTR is correctly compiled by running `./superstr --help`; you should see:

```
Usage: superstr [-s] [--help] [--version] [--mode=<fastq|bam>] [-o myfile] <file> [<file>]...
Rapid STR characterisation in NGS data.
  --help                    display this help and exit
  --version                 display version info and exit
  --mode=<fastq|bam>        type of input data
  -s, --stream              run on named streams, not files (see manual for instructions)
  -o myfile                 output directory
  <file>                    input files (or names of pipes in stream mode)
```

If you run into issues with installation, check the HTSLIB_ROOT variable, and then try deleting the cmake cache (at superSTR/CMakeCache.txt) and re-running the install process from step 5.

### post-processing

The easiest way to install superSTR's python dependencies is via the conda package manager - simply run `conda env create -f environment.yml` while working in the superSTR/ root.

Failing that, the current list of dependencies is:

> - python =3.8.3
>  - numpy =1.20.1
>  - matplotlib =3.4.1
>  - pandas =1.2.2
>  - bashtage::arch =4.15
>  - mpmath =1.2.1
>  - scipy =1.6.1
>  - statsmodels =0.12.2
>  - tqdm =4.58.0
>  - seaborn =0.11.1

## Installation via Docker

Docker and Singularity files in this repository are provided as examples and have not been exhaustively tested. 

All care has been taken in their preparation, but no warranty is made as to their security; you use them at your own risk; we strongly recommend using a tool like [Snyk](snyk.io) to check the security of this image prior to building and running it, or getting your local IT team to check. The Dockerfile in this repository is kept up to date as Snyk reports security fixes.

These files currently contain only the components required for read processing steps; you will need to install the postprocessing code as outlined above.

Docker: [superSTR/Dockerfile](Dockerfile)
