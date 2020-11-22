# Outlier detection

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

## Metrics

1. Library size-normalised count at max read length
2. superSTR summary score

## Outlier detection methods

**ABC-bootstrapped 95th quantile estimator**
**Count vector autoencoder**
