# superSTR

23/11: superSTR is currently being loaded into this repo for an initial public release; there will likely be some minor changes to code and documentation over the next 48 hours.

A lightweight, alignment-free utility for detecting repeat-containing reads in short-read WGS, WES and RNA-seq data. A preprint describing superSTR is in preparation.

The C99 component of superSTR is a modified and extended version of [mreps](https://github.com/gregorykucherov/mreps) by Roman Kolpakov, Ghizlane Bana and Gregory Kucherov. Full details of mreps can be found at http://mreps.univ-mlv.fr/ and in its accompanying paper; R. Kolpakov, G. Bana, and G. Kucherov, mreps: efficient and flexible detection of tandem repeats in DNA, Nucleic Acid Research, 31 (13), July 1 2003, pp 3672-3678.

Full details of libraries used in superSTR can be found in [the Acknowledgements file](docs/ACKNOWLEDGEMENTS.md).

## Basic superSTR operations

This section describes how to run basic superSTR analysis with a minimum of fuss on human genomic samples.

[Installing superSTR](docs/INSTALL.md)

1) [Processing FASTQs and BAM files](docs/PROCESSING.md)
2) [Post-processing of sets of samples](docs/POSTPROC.md)
3) [Outlier detection](docs/OUTLIERS.md)
4) [Visualisation](docs/VISUALISATION.md)

## Extending and advanced superSTR

This section contains information necessary to run superSTR outside of human genomic samples.

* [Compression ratio thresholds and how to set them](docs/THRESHOLD.md)
