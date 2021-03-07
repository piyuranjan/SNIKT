# SNIKT

Identify adapter/systemic contamination in metagenomic sequencing DNA/RNA data.

### Installation

SNIKT is only available on Unix-like platforms (Linux, Mac OS and WSL). The easiest way to use SNIKT, in our opinion, is to use the Conda package management for set up. However, being a wrapper written in R, this program can also be used with a local R installation given it has access to `seqtk` and some R libraries. Both of these setup options are explained below.

#### Option 1: Installation via Bioconda

For this setup process to work, you will need to install the Conda package manager and setup Bioconda repositories. Instructions for setting up Bioconda can be found on the [Bioconda installation page](https://bioconda.github.io/user/install.html). Only setup steps 1 and 2 are needed for this. Once you have the `conda` set up and working, please start a new session.

We are considering making this software available as a conda recipe in the future which will simplify the installation to almost a single command. Until then, please follow these steps to set up SNIKT via `conda`.

Set up a new environment and install dependencies. You can make your own choice but here we explain these steps with `envSnikt` as our environment name.
```
$ conda create -n envSnikt r-tidyverse r-gridExtra r-docopt r-lubridate seqtk
```
Download the SNIKT code from the GitHub repository and place it under the environment we just created. For this purpose, we are assuming that we had set up `miniconda` in the location `~/miniconda3/`. Please make sure to use the right path to your `conda` installation. This path shows up during `bioconda` setup as well. After keeping the `snikt.R` code in the environment, we provide it execute permissions.
```
$ wget https://github.com/piyuranjan/SNIKT/raw/main/snikt.R -O ~/miniconda3/envs/envSnikt/bin/snikt.R
$ chmod +x ~/miniconda3/envs/envSnikt/bin/snikt.R
```
And that's it for the setup with this option. Activate the environment and get the power of Wolverine's claws in your hands!
```
$ conda activate envSnikt
$ snikt.R --version
SNIKT 0.3.0
```
If you ever need to uninstall, you can directly remove the environment with the following.
```
conda remove -n envSnikt --all
```

#### Option 2: Using a local R installation

While the `conda` installation is easier, we recognize that many users may want to use their local R installation instead of pulling a fresh copy of R and its libraries as a different environment. So, to setup SNIKT in this way, you will need to first get the program `seqtk` from the [Seqtk GitHub page](https://github.com/lh3/seqtk). A slightly older version of Seqtk is also available in the `apt` package manager for Ubuntu Linux.

Next, you can install the packages in your R (we recommend R >= 4.0) library. To do this, you can do the following in your R terminal.
```
$ R
R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
...

> install.packages( c("tidyverse","gridExtra","docopt","lubridate"), dependencies = TRUE)
```

Now, that the dependencies are set, we can download SNIKT.
```
$ git clone https://github.com/piyuranjan/SNIKT.git
```
After this download, you can either put the downloaded folder on your system `PATH` or copy the `snikt.R` code to a location that is already on your `PATH`.


### Dependencies

This program has been tested with the following dependencies and their versions.

- seqtk >= 1.3 [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)
- R >= 4.0.3 [https://www.r-project.org/](https://www.r-project.org/) [direct download](https://cloud.r-project.org/) with the following libraries.
  - tidyverse 1.3.0
  - grid 
  - gridExtra 2.3
  - docopt 0.7.1
  - lubridate 1.7.10

### Full command line help and usage

```
SNIKT: FastQ QC and sequence over-representation check.
       A wrapper around seqtk to plot per-position nucleotide composition
       for finding and trimming adapter contamination in fastq reads.
       Also filters reads by a length threshold.
Authors: Piyush Ranjan, Christopher Brown

For first-time users, interactive mode is recommended.
For detailed help and examples, please visit
https://github.com/piyuranjan/SNIKT

Location: scriptPath

Usage:
  snikt.R [options] [--] <fastq>
  snikt.R <fastq>  # Interactive # Easiest
  snikt.R [--zoom5=<nt> --zoom3=<nt>] <fastq>  # Interactive
  snikt.R [(--trim5=<nt> --trim3=<nt>) | --notrim] <fastq>
  snikt.R [--illumina] [-n] <fastq>

Input:
  <fastq>               Sequence file in fastQ format with exts: .fq, .fq.gz,
                          .fastq, .fastq.gz

Options:
  Presets:
  --illumina            This presets options that are better for short-read
                          Illumina datasets.
                          Sets: -f 0 -Z 50 -z 50
                          Defaults are configured for long-read Nanopore fastq.

  Graphing:
  --hide=<frac>         Hide the composition tail by a fraction of total bases.
                          Significantly improves speed, removes end-tail (3')
                          distortion for variable length read sets.
                          [range: 0..1] [default: 0.01]
  -s, --skim=<num>      Use top num reads for pre- or no-trim graphs. This
                          improves speed. No effect on post-trim graphs.
                          Use 0 to disable skimming and utilize all reads.
                          [range: 0..maxFastqReads] [default: 10000]
  -Z, --zoom5=<nt>      Zoom-in from aligned 5' beginning to nt bases.
                          [range: 1..maxSeqLen] [default:300]
  -z, --zoom3=<nt>      Zoom-in from aligned 3' ending to nt bases.
                          [range: 1..maxSeqLen] [default:100]
  QC:
  -n, --notrim          Disable positional trimming; useful for short-read data
                          Takes precedence over and sets: -T 0 -t 0
  -T, --trim5=<nt>      Trim nt bases from aligned 5' side.
                          [range: 0..(maxSeqLen-trim3)] [default: interactive]
  -t, --trim3=<nt>      Trim nt bases from aligned 3' side.
                          [range: 0..(maxSeqLen-trim5)] [default: interactive]
  -f, --filter=<nt>     Filter (drop) reads with length < nt after any trimming.
                          [range: 0..maxSeqLen] [default:500]
  IO:
  -o, --out=<prefix>    Prefix for output files [default: fastqNoExtension]
  -w, --workdir=<path>  Path to generate QC file and report. [default: ./]
  --gzip                If fastq file is gzipped. Autodetected normally
                          using the file extension. For large datasets, prior
                          decompression of fastq may be faster. Only gzip is
                          supported within; prior decompression needed for any
                          other method.
  -k, --keep            Keep intermediate (temporary) directory.
  
  Generic:
  -h, --help            Show help and exit 0.
  -v, --verbose         Enable status messages.
  --debug               Debug with traceback; enables -v.
  --version             Show version and exit 0.
```

### Citation

A manuscript is under preparation for this project. Until then, please feel free to cite this project with its GitHub URL.
