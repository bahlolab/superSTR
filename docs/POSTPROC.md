# Post-processing superSTR C output

superSTR generates a single file (suffixed per_read.txt.gz) for each sample that contains one line per repeat-containing read. 

The next step in superSTR analysis is to produce a set of per-sample and per-motif files, which is achieved using the Python multiparse.py script.


## Input

The post-processing step requires superSTR per-sample files, and a manifest file that contains label and path information.

**Manifest file**

Manifest files are derived from [ExpansionHunter deNovo manifests](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/06_Merging_profiles.md); they are a four-column file with the following information:

| sample ID | case/control | file path | group label |

The first two columns can be used in an EHDN manifest file.

## Running post-processing

**Requirements:** multiparse.py was developed and tested with Python >= 3.7. It is untested in versions of Python < 3.7. It requires at least two processes (one or more reader processes and a writer process). Your operating system may be able to schedule these if you have less than two CPUs (*e.g.* if you are running on a small 1 CPU cloud instance), but this is not currently tested.

**With a manifest file (4 processes, 150nt reads):** python multiparse.py -m manifest.tsv --output output_dir/ -@ 4 -r 150

## Post-processing output

Post-processing produces a directory structure with the following structure:

```
analysis/
└───samples/
│   │   per_read.txt.gz   
└───motifs/
│   │   3mers/
|   │   │   AAC.csv
|   │   │   ...
│   │   4mers/
|   │   │   AAAC.csv
|   │   │   ...
|...
```
