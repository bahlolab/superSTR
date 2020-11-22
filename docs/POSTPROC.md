# Post-processing superSTR C output

superSTR generates a single file (suffixed per_read.txt.gz) for each sample that contains one line per repeat-containing read. 

The next step in superSTR analysis is to produce a set of per-sample and per-motif files, which is achieved using the Python multiparse.py script.


## Input

The post-processing step requires superSTR per-sample files, and one of two input structures; either a manifest file that contains label and path information, or a directory path.

**Manifest file**

Manifest files are based on [ExpansionHunter deNovo manifests](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/06_Merging_profiles.md); they are a four-column file with the following information:

| sample ID | case/control | file path | group label |

Deleting the last column with your preferred editor or utility will produce an EHdN manifest file.

**Input directory**

If a directory is provided, superSTR expects it to have the following format (the default produced by the C superSTR method when provided with `output/sample<num>/` as the output folder).

```
input/
└───sample1/
│   │   per_read.txt.gz   
└───sample2/
│   │   per_read.txt.gz  
|...
```

## Running post-processing

**Requirements:** multiparse.py was developed and tested with Python 3.7; it will not work on Python 2, and is untested in versions of Python < 3.7. 

