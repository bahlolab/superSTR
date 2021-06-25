# Motif screening

Motif screening computes the Mann-Whitney statistic to compare the distribution of information scores between cases and controls for each motif. 

The exact p-value for the statistic is estimated using permutation testing.

## Time and resource requirements

Permutation testing is highly resource-intensive, and the run-time for this analysis depends on how many motifs have been detected and how many conditions are being tested. 

The analysis of 3-6mers in the SCA3 RNA-seq data (12 cases, 12 controls) described in the superSTR manuscript and from which data has been taken for the tutorial completed in <3 minutes on a 2015 Macbook Pro. Analysis of the large number of different ICD codes in the 50,000 samples in the UK Biobank required significant parallelisation and HPC resources.

# Basic usage

`python superSTR/Python/screen.py -m manifest.tsv -i input_path/ -o output.tsv --controllab Control --min_thresh N`

The `input_path/` folder should point to a directory with the following structure, as produced by the postprocessing step:

```
input_path/
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

The `min_thresh` parameter is used to control the number of samples that must be in a group for that group to be tested - for example, a min-thresh of 1 will attempt to compare case labels with 1 sample to the control set. It is important to note that the Mann-Whitney statistic compares distributions; it is more-appropriate to run outlier detection should there be a small number of samples for a particular condition.
