# MetagenomeAdapterIdentify
Identify adapter/systemic contamination in metagenomic sequencing DNA/RNA data.

### Dependencies
seqtk >= 1.3 [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

R >= 3.5.1 [https://www.r-project.org/](https://www.r-project.org/) [direct download](https://cloud.r-project.org/)  

Libraries:
- tidyverse
- grid
- gridExtra
- argparser
- lubridate

### Command line help and usage
```
readAdapterIdentify.R [--help] [--verbose] [--filter FILTER] [--bzoom BZOOM] [--ezoom EZOOM] [--outfile OUTFILE] fastq

FastQ QC and sequence over-representation check.
Calculate and plot per-position nucleotide composition, for finding
adapter contamination/over-representation in sequences.
Authors: Piyush Ranjan, Christopher Brown

For more info, please check: https://github.com/piyuranjan/MetagenomeAdapterIdentify

positional arguments:
  fastq                     FastQ file for composition analysis; required

flags:
  -h, --help                show this help message and exit
  -v, --verbose             Enable status messages [default: FALSE]

optional arguments:
  -f, --filter FILTER       [0..1] Filter the composition tail by a fraction [default: 0.005]
  -b, --bzoom BZOOM         [1..maxSeqLen] NT to zoom in from aligned 5' begin [default: 300]
  -e, --ezoom EZOOM         [1..maxSeqLen] NT to zoom in from aligned 3' end [default: 100]
  -o, --outfile OUTFILE     [file.ext] File for saving graphs [default: fastqNoExt-fastqCompositions.png]
```
