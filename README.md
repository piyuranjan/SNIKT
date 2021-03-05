# SNIKT

Identify adapter/systemic contamination in metagenomic sequencing DNA/RNA data.

### Dependencies

seqtk >= 1.3 [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

R >= 4.0.3 [https://www.r-project.org/](https://www.r-project.org/) [direct download](https://cloud.r-project.org/)  

Libraries:
- tidyverse
- grid
- gridExtra
- docopt
- lubridate

### Command line help and usage

```
SNIKT: FastQ QC and sequence over-representation check.
       A wrapper around seqtk to plot per-position nucleotide composition
       for finding, length trimming and filtering fastq sequences.
       Removes sequence-end and adapter contamination.
Authors: Piyush Ranjan, Christopher Brown

For detailed help and examples, please visit
https://github.com/piyuranjan/SNIKT

Location: ./snikt.R

Usage:
  snikt.R [options] [--] <fastq>
  snikt.R <fastq>  # Interactive
  snikt.R [--zoom5=<nt> --zoom3=<nt>] <fastq>  # Interactive
  snikt.R [(--trim5=<nt> --trim3=<nt>) | --notrim] <fastq>
  snikt.R [--illumina] [-n] <fastq>

Input:
  <fastq>               Sequence file in fastQ format with ext .f[ast]q[.gz]

Options:
  Presets:
  --illumina            This presets options for short-read Illumina datasets.
                          Sets: -f 0 -Z 50 -z 50
                          Defaults are set for long-read Nanopore sequences.

  Graphing:
  --hide=<frac>         Hide the composition tail by a fraction of total bases.
                          Significantly improves speed, removes end-tail (3')
                          distortion for variable length read sets.
                          [range: 0..1] [default: 0.01]
  -s, --skim=<num>      Skim fastq with top num reads for pre- or no-trim
                          graphs. Improves speed. No effect on post-trim graphs.
                          Use 0 to disable skimming and utilize all reads.
                          [range: 0..maxFastqReads] [default: 10000]
  -Z, --zoom5=<nt>      NT to zoom-in from aligned 5' beginning
                          [range: 1..maxSeqLen] [default:300]
  -z, --zoom3=<nt>      NT to zoom-in from aligned 3' ending
                          [range: 1..maxSeqLen] [default:100]
  QC:
  -n, --notrim          Disable positional trimming; useful for short-read data
                          Takes precedence over and sets: -T 0 -t 0
  -T, --trim5=<nt>      NT to trim from aligned 5' side
                          [range: 0..(maxSeqLen-trim3)] [default: interactive]
  -t, --trim3=<nt>      NT to trim from aligned 3' side
                          [range: 0..(maxSeqLen-trim5)] [default: interactive]
  -f, --filter=<nt>     Filter (drop) reads with length < nt after any trimming.
                          [range: 0..maxSeqLen] [default:500]
  IO:
  -o, --out=<prefix>    Prefix for output files [default: fastqNoExt]
  -w, --workdir=<path>  Path to generate QC file, report. [default: ./]
  --gzip                If fastq file is gzipped. Autodetected normally
                          using the file extension. For large datasets, prior
                          decompression of fastq may be faster.
  -k, --keep            Keep intermediate (temporary) directory.

  Generic:
  -p, --proc=<threads>  Num proc threads for parallel trimming. All other
                          actions are still serial. Requires GNU parallel.
                          [range: 1..maxCores] [default: 1]
  -h, --help            Show help and exit 0.
  -v, --verbose         Enable status messages.
  --debug               Debug with traceback; enables -v.
  --version             Show version and exit 0.
```
