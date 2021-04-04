# Acknowledgements

superSTR uses code from a variety of sources as libraries or by forking and modifying projects. We are grateful to the developers and maintainers of these packages; this list attempts to capture the dependencies used in superSTR in addition to citations elsewhere.

# C

## mreps - C

The C99 component of superSTR is a modified and extended version of mreps by Roman Kolpakov, Ghizlane Bana and Gregory Kucherov. Full details of mreps can be found at http://mreps.univ-mlv.fr/ and in its accompanying paper; R. Kolpakov, G. Bana, and G. Kucherov, mreps: efficient and flexible detection of tandem repeats in DNA, Nucleic Acid Research, 31 (13), July 1 2003, pp 3672-3678.

## argtable - C

superSTR's argument parsing uses the [argtable](https://github.com/argtable/argtable3) library, and its argument parsing is based on the argtable examples/tutorials.

## queue.h - C

superSTR contains a copy of the [BSD queue library](https://www.freebsd.org/cgi/man.cgi?query=queue&sektion=3&manpath=FreeBSD+12.2-RELEASE+and+Ports). License information for BSD queue can be found in the header of queue.h.

# Python

## Autoregressive Conditional Heteroskedasticity (ARCH)

superSTR's outlier detection method uses the IIDBootstrap method from the [ARCH library](https://pypi.org/project/arch/). Details of this package, including license information, are available [here](https://github.com/bashtage/arch).

## mpmath 

superSTR uses the [mpmath library](https://github.com/fredrik-johansson/mpmath) for arbitrary precision mathematics in calculating exact p-values during permutation testing, particularly for accurate Gauss-Legendre quadrature. A full list of contributors to mpmath are listed in the project's Github page.

