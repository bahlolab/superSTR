# Outlier detection

23/11: Note: the current outlier.py script only implements the ABC-bootstrapped 95th quantile estimator due to a bug in the install process for the libraries necessary for the other methods. This is being worked on. 

Outlier detection within superSTR output is context dependent, but simple implementations of three methods and two built-in metrics are provided in this library along with instructions for extending these methods or implementing your own outlier detection methods.

## Requirements

The provided outlier detection code requires installation of some python libraries, which I recommend you install using conda (or venv).

```
conda create --name superSTR python=3.7
pip install -U pip
pip install arch
...
conda activate superSTR
```

## Time and resource requirements

Outlier detection methods may be highly resource-intensive, or produce outputs for very large numbers of motifs. To help with managing this load we provide a set of (optional) flags that allow you to limit the outlier detection process to certain motifs by identity or length. 

They are:

`--pathogenic`

Limits outlier detection to known-pathogenic motifs sourced from Bahlo, M. *et al*, [Recent advances in the detection of repeat expansions with short-read next-generation sequencing](https://f1000research.com/articles/7-736/v1), F1000 Research (2018).

`--max_motif n`

Limits outlier detection to motifs with a motif length <= n.


## Metrics

1. Library size-normalised count at max read length
2. superSTR summary score

## Outlier detection methods

**ABC-bootstrapped 95th quantile estimator**

