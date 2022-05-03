# superSTR

A lightweight, alignment-free utility for detecting repeat-containing reads in short-read WGS, WES and RNA-seq data.

The C99 component of superSTR uses and incorporates a modified version of [mreps](https://github.com/gregorykucherov/mreps) by Roman Kolpakov, Ghizlane Bana and Gregory Kucherov. Full details of mreps can be found at http://mreps.univ-mlv.fr/ and in its accompanying paper; R. Kolpakov, G. Bana, and G. Kucherov, mreps: efficient and flexible detection of tandem repeats in DNA, Nucleic Acid Research, 31 (13), July 1 2003, pp 3672-3678.

Full details of libraries used in superSTR can be found in [the Acknowledgements file](docs/ACKNOWLEDGEMENTS.md).

## Basic superSTR operations

This section describes how to run basic superSTR analysis with a minimum of fuss on human genomic samples.

[Installing superSTR](docs/INSTALL.md)

1) [Processing FASTQs and BAM files](docs/PROCESSING.md)
2) [Post-processing of sets of samples](docs/POSTPROC.md)
3) [Outlier detection](docs/OUTLIERS.md)
4) [Motif screening](docs/SCREENING.md)
5) [Visualisation](docs/VISUALISATION.md)

We also provide a detailed [example RNA-seq analysis](docs/EXAMPLE.md) based on the SCA3 data used in the superSTR manuscript to illustrate an end-to-end superSTR analysis.

## Upcoming changes:

* A few modifications to improve UI in the code will be implemented to standardise flags across utilities.
* mreps' resolution parameter will be exposed to the user in the next update.

# Licensing

superSTR is released under the GNU General Public License, v2.

# Acknowledgements:

* @jbedo for bugfixes and helping streamline the Cmake files used to compile superSTR.
* @christopher-schroeder for developing conda installation methods for superSTR and bugfixes.
