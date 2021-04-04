# Installation:

## Requirements:

superSTR has not been tested on (and does not currently support) Microsoft Windows; [a containerised version](#installation-via-docker) of superSTR is available but is untested on that platform. 

The following instructions have been tested on:

* Debian 10.8 "buster"
* Ubuntu 20.04 LTS "Focal Fossa"
* Centos 8
* Apple OSX 10.14.6 "Mojave" 

Please make sure the following are installed and available on your system prior to installing superSTR:

* [htslib 1.9](https://github.com/samtools/htslib)
* [zlib](https://zlib.net/)
* [Cmake](https://cmake.org/install/) (>3.10)
* [curl](https://curl.se/download.html) and [libcurl](https://curl.se/libcurl/)
* [gcc](https://gcc.gnu.org/install/)
* [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

## Installation steps - read processing:

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

## Installation steps - postprocessing

We recommend the use of the conda package manager. 

To install all dependencies for superSTR using conda (or miniconda), use `conda env create -f superSTR.yml` from the repository main directory. This will create a superstr environment with all dependencies installed.

## Installation via Docker

Docker and Singularity files in this repository are provided as examples. All care has been taken in their preparation, but no warranty is made as to their security; you use them at your own risk; we strongly recommend using a tool like [Snyk](snyk.io) to check the security of this image prior to running it, or getting your local IT team to check.

These files currently contain only the components required for read processing steps; you will need to install the postprocessing code as outlined above.

Docker: [superSTR/Dockerfile](Dockerfile)
