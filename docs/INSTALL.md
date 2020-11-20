# Installation:

Requirements:

* htslib 1.9
* Cmake >3.10
* zlib/libzlib
* curl/libcurl

1) Clone this repo using `git clone https://github.com/bahlolab/superSTR`.
2) Change into the superSTR directory with `cd superSTR`
3) Set the the HTSLIB_ROOT environment variable with `export HTSLIB_ROOT=<path_to_your_htslib_installation>`
4) Change into the C source directory with `cd C`
5) Run `cmake . && make .`
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

# Installation - Docker

#TODO: Docker
