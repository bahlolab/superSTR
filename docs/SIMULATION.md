# Simulation and compression threshold setting

There are two simulation methods that can be used for determining appropriate thresholds for use in superSTR.

* simulation.py (documented here) creates simulated, pseudorandom reads and evaluates threshold performance under different measures.
* ART-based simulated read analysis

## simulation.py

simulation.py operates in two modes which sare certain parameters:

### Common Parameters


*Output directories - required, all modes*

`--output`: path to a directory to contain output - eg `--output ./`

*Distribution flags - all modes*

`--hsa`: Homo sapiens background CD proportion of 40.87% per https://doi.org/10.1186/s13104-019-4137-z

`HSA_DNA_PROBABILITIES = [0.20435, 0.20435, 0.29565, 0.29565]`

`--pfa`: AT content of intergenic regions of P falciparum per doi:10.1093/nar/gkw1259

PFA_DNA_PROBABILITIES = [0.068, 0.068, 0.432, 0.432]`

`--sco`: approximation of the S. coelicolor C/G:A/T ratio per https://doi.org/10.1038/417141a

`SCO_DNA_PROBABILITIES = [0.3606, 0.3606, 0.1394, 0.1394]`

`--even`: even probabilities of C, G, T and A.

`EVEN_DNA_PROBABILITIES = [0.25, 0.25, 0.25, 0.25]`

A convenience option `--all`, is provided, equivalent to `--hsa --pfa --sco --even`

## Simulation mode: --simulation

Randomly samples repeat reads of various types (both partial and pure, along with 'background' sequence), then plots distributions of the compression ratios.

Requires specification of a threshold for use with plotting (`-t` option).

The type of background sequence may be selected by the user with one of the following options:

Use a pseudorandom sequence from a predefined nucleotide distribution:

`simulation.py --simulation -t <threshold> --output <output_path> <distribution flag(s)>`

Provide a custom reference as a gzipped FASTA file and run with:

`simulation.py --simulation -t <threshold> --reference_path <reference_path> --output <output_path>`

Download sequence from ENSEMBL GRCh37, using contig lengths provided by UCSC:

`simulation.py --simulation -t <threshold> --ensembl --output <output_path>`

## Performance

Performance mode randomly samples repeat reads of various types (both partial and pure repeat, along with pseudorandom sequence), then prints out the optimal threshold under the ROC and precision+recall metrics under different predefined GC/nucleotide distributions.

For default read lengths:

`simulation.py --performance --all --output <output_path>`

For a custom read length:

`simulation.py --performance --all --output <output_path> -r <read_length>`