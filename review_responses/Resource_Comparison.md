# Comparison of resource usage
Created: 2022-04-08  
Author: Piyush Ranjan

---

[GitHub Location](https://github.com/piyuranjan/SNIKT/blob/main/review_responses/Resource_Comparison.md): `piyuranjan/SNIKT/tree/main/review_responses/Resource_Comparison.md`
Local Location: `$PWD/Resource_Comparison/` (disregard on GitHub)

In response to reviewer comments, we compared the performance of SNIKT with three popular tools for sequencing read quality control (QC) - [NanoQC](https://github.com/wdecoster/nanoQC), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [NanoFilt](https://github.com/wdecoster/nanofilt). NanoQC and NanoFilt are from the [NanoPack](https://github.com/wdecoster/nanopack) ([DOI](https://doi.org/10.1093/bioinformatics/bty149)) suite of tools that are commonly used for Nanopore long-read QC. FastQC is a QC program very commonly used for Illumina and PacBio reads, however, it can be modified to work with Nanopore long reads with varying degree of success.

In this analysis, we draw comparisons over resource utilization for the different tools and their effectiveness in QC for ONT long-read data. We categorize them as the following.
- Comparison of QC reporting: NanoQC, FastQC and SNIKT
- Comparison of trimming: NanoFilt and SNIKT

<br>

## Test setup

### Dataset
The dataset used for this comparison is the same SRA run ([Garg et al. 2021](https://doi.org/10.1128/MRA.00299-21)) that is used in the main paper. This sequencing run is a WGS long-read Nanopore run with 9.4.1 pore chemistry for HMW DNA extraction from *Candida albicans* CHN1 which was prepared with the SQK-RAD004 rapid DNA library preparation. Following are basic summary statistics from this fastq using a custom script.
```_
$ summarizeFastq.pl SRR13441294.fastq.gz
#Dataset        #Reads  #Bases  #MinLen #MaxLen #AvgLen #AvgQ   #AvgA   #AvgC   #AvgG   #AvgT   #AvgN
SRR13441294.fastq.gz    412984  4900224854      1       166478  11865.41        21.3    32.7    16.9    17.2    33.2    0.0
```
| Dataset              | Reads  | Bases      | MinLen | MaxLen | AvgLen   | AvgQ | AvgA | AvgC | AvgG | AvgT | AvgN |
| -------------------- | ------ | ---------- | ------ | ------ | -------- | ---- | ---- | ---- | ---- | ---- | ---- |
| SRR13441294.fastq.gz | 412984 | 4900224854 | 1      | 166478 | 11865.41 | 21.3 | 32.7 | 16.9 | 17.2 | 33.2 | 0.0  |

The following adapter sequences known to be used for the [Rapid Sequencing kit family](https://community.nanoporetech.com/technical_documents/chemistry-technical-document/v/chtd_500_v1_revae_07jul2016/rapid-sequencing-kit-family) (SQK-RAD004) were gathered from the ONT community website.
Adapter Y top, 61 nt
```_
5'- GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT -3'
```
Adapter Y bottom, 22 nt
```_
5'- GCAATACGTAACTGAACGAAGT -3'
```
These sequences were kept in fasta and tab-separated formats. However, according to the ONT sequencing chemistry, only Adapter_Y_top is supposed to be seen in this read set.

### Device
All tests were done on a 12 core Windows laptop running WSL2 with Ubuntu 20.04 as the Linux distribution.
```_
$ lscpu | head -n 13
Architecture:                    x86_64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
Address sizes:                   39 bits physical, 48 bits virtual
CPU(s):                          12
On-line CPU(s) list:             0-11
Thread(s) per core:              2
Core(s) per socket:              6
Socket(s):                       1
Vendor ID:                       GenuineIntel
CPU family:                      6
Model:                           158
Model name:                      Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz

$ uname -a
Linux MSL00332 5.10.16.3-microsoft-standard-WSL2 #1 SMP Fri Apr 2 22:23:49 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux
```

### NanoQC and NanoFilt installation
Both NanoQC and NanoFilt were installed in a single Conda environment.
```_
$ mamba create -n nanopack nanoqc nanofilt
$ mamba activate nanopack
$ nanoQC -v
NanoQC 0.9.4
$ NanoFilt -v
NanoFilt 2.8.0
```

### FastQC installation and setup
FastQC was setup with modified memory buffer, so it could work with moderately long ONT reads.
FastQC v0.11.9 was setup in a Conda environment.
```_
$ mamba create -n fastqc fastqc
$ mamba activate fastqc
$ fastqc --version
FastQC v0.11.9
```
FastQC requires modification of the memory buffer to allow for larger heap space for Java. Until that is done, FastQC fails to analyze most long-read sets with the following error. Even after this modification, ultra-long read sets still cause memory errors for FastQC.
```_
$ fastqc SRR13441294.fastq
Started analysis of SRR13441294.fastq
Exception in thread "Thread-1" java.lang.OutOfMemoryError: Java heap space
        at uk.ac.babraham.FastQC.Utilities.QualityCount.<init>(QualityCount.java:33)
        at uk.ac.babraham.FastQC.Modules.PerBaseQualityScores.processSequence(PerBaseQualityScores.java:141)
        at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:89)
        at java.base/java.lang.Thread.run(Thread.java:834)
```
Java memory buffer modification was performed in `~/.miniconda3/envs/fastqc/opt/fastqc-0.11.9/fastqc` as below.
```_
$ cd .miniconda3/envs/fastqc/bin
$ ll fastqc
lrwxrwxrwx 1 pr pr 27 Apr 13 12:01 fastqc -> ../opt/fastqc-0.11.9/fastqc*
$ cd ../opt/fastqc-0.11.9/
$ nano fastqc
```
Following code was modified changing single-thread instance memory from 250 MB to 10 GB.
```perl
if ($threads) {
	if ($threads < 1) {
		die "Number of threads must be a positive integer";
	}
	
	push @java_args ,"-Dfastqc.threads=$threads";
	my $memory = 250 * $threads;
	my $stack = 200 * $threads;
	unshift @java_args,"-Xmx${memory}m";
	unshift @java_args,"-Xms${memory}m";
}
else {
	# unshift @java_args,'-Xmx250m';
	unshift @java_args,'-Xmx10g';
	# unshift @java_args,'-Xms200m';
	unshift @java_args,'-Xms10g';
}
```
In testing the installation, these parameters seemed to work for the fastq dataset used in this analysis.

### SNIKT installation
SNIKT was installed via Conda. However, at the time of this writing version `v0.4.3` was not released on the Anaconda. So, a local copy of the script was created and used, while all dependencies came from the Conda environment.
```_
$ mamba create -n snikt snikt
$ mamba activate snikt
$ snikt.R --version
SNIKT 0.4.2
$ wget https://raw.githubusercontent.com/piyuranjan/SNIKT/main/snikt.R
$ chmod 755 snikt.R
$ ./snikt.R --version
SNIKT 0.4.3
```

### Timing of commands
Timing of the processes were done with `/usr/bin/time`. A specific format of the command was aliased as `recordStats` for repeated and consistent usage.
```_
$ alias recordStats
alias recordStats='/usr/bin/time -f "\nCommand stats:\nProc:\tElapsed Time = %E,\tAvg CPU = %P,\nMem:\tAvg Total Mem = %KKB,\tPeak Mem = %MKB,\nDisk:\tIn = %I,\tOut =%O\nExit Status: %x"'
```

<br>

## Comparison of QC reporting

### FastQC with `--adapter` execution
FastQC was run with the supplemented adapters to prepare a QC report.
```_
$ recordStats fastqc -a ../adapters.txt SRR13441294.fastq.gz
Started analysis of SRR13441294.fastq.gz
Approx 5% complete for SRR13441294.fastq.gz
...
Approx 95% complete for SRR13441294.fastq.gz
Analysis complete for SRR13441294.fastq.gz

Command stats:
Proc:   Elapsed Time = 3:53.67, Avg CPU = 101%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 7286240KB,
Disk:   In = 9885208,   Out =4248
Exit Status: 0
```
Following report was generated. Download locally to see.  
[FastQC Report: ./\_attach/SRR13441294_fastqc.html](_attach/SRR13441294_fastqc.html)

#### FastQC result summary
- While this process finished in under 4 minutes on a single thread, it took 7.3 GB RAM memory. The documentation suggests that this memory utilization is per thread, which means higher thread count could easily push this requirement to unsustainable values. This program has been seen to multi-thread very efficiently on short-read Illumina datasets with a moderate RAM impact, however, that observation does not apply for long-read datasets.
- While a custom file for adapters was used, the adapter detection in the QC report did not seem to work. This could be due to high error rate in ONT data which often presents difficulty in adapter alignment. This is especially true for packages like FastQC which are optimized based on accuracies seen on Illumina instruments.

### NanoQC default `--minlen` execution
NanoQC was executed without `--minlen` as the description seems to suggest it will filter reads before making a QC report.
```_
$ recordStats nanoQC SRR13441294.fastq.gz
BokehDeprecationWarning: plot_width and plot_height was deprecated in Bokeh 2.4.0 and will be removed, use width or height instead.

Command stats:
Proc:   Elapsed Time = 18:31.92,        Avg CPU = 100%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 1759732KB,
Disk:   In = 0, Out =144
Exit Status: 0
```
Following report was generated. Download locally to see.  
[NanoQC Report: ./\_attach/nanoQC.html](_attach/nanoQC.html)

#### NanoQC result summary
- The default for `--minlen` was taken as 200 and the method did filter reads smaller than that. It made graphs for a 100 bp window from each end. This means that this package necessitates a coupling between the viewing window and the read length filter removing a true representation of the original file before plotting it. This could be detrimental if the fastq file has a higher fraction of small reads. The log suggested that it used 411060 reads instead of the original 412984 which leaves 1924 reads outside of the initial representation.
- The graphing step took about 18.5 minutes, using 1.7 GB RAM memory while working over a single thread. This application is not multi threaded (neither is SNIKT).
- While addressable, some graphs are not labelled.
	- While the graphs are interactive, the data drawn on them is fixed. It limits the interactivity to only zooming in.
	- The frequency length graph is set as a separate graph in NanoQC report. In SNIKT, this graph is embedded with compositions and Phred score on the same x-axis, which allows the user to make a simultaneous decision on how much throughput they would waste at different levels of trimming along with choosing an appropriate criteria using the nucleotide compositions and the Phred scores.

### SNIKT `--notrim` execution
SNIKT was run with its `-n, --notrim` parameter so it makes a report for initial analysis of the adapter and low quality bases in the fastq. This behavior is analogous to how NanoQC works.  
While SNIKT defaults to using top 10,000 reads for its graphs for a significantly faster assessment, to make the comparison fair, we also used the `-s, --skim`  parameter which enables the user to force SNIKT to use the entire file for graphing.  
Default parameters for 5' and 3' zoom criteria were used but they can be modified to the values more suitable for adapter assessment using `-Z, --zoom5` and `-z, --zoom3` flags. Filtering reads at this step is by design disabled, so using `-f, --filter` will have no effect as SNIKT is only making an estimation and we think it is important to be able to get an accurate representation of the fastq.
```_
$ recordStats ./snikt.R --skim=0 --notrim --out=SRR13441294_notrim SRR13441294.fastq.gz
SNIKT contamination report is available in: /home/pr/WD/WD2/SRR13441294_notrim.html

Command stats:
Proc:   Elapsed Time = 4:18.89, Avg CPU = 132%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 757896KB,
Disk:   In = 160488,    Out =95088
Exit Status: 0
```
Following report was generated. Download locally to see.  
[SNIKT pre-trim report: ./\_attach/SRR13441294_notrim.html](_attach/SRR13441294_notrim.html)

#### SNIKT result summary
- Unlike NanoQC, SNIKT took the entire fastq into consideration without any filters. The default x-axis window sizes are optimized for long-read shotgun sequencing. They work very well in showing the poor quality bases and adapter contamination in an easy view, so that the user can determine where to make their trims. The left and the right panels also demonstrate negligible loss in number of reads (using read frequency) when making trims at, what we see appropriate, 100 bp from 5' and 20 bp from 3' for this set. The graph is an image that can be opened in a bigger view, if desired, in a new tab or can be saved to zoom in.
- The graphing step took about 4.3 minutes (NanoQC - 18.5 min, FastQC - 4 min) to complete, using 760 MB in RAM memory (NanoQC - 1.7 GB, FastQC - 7.3 GB) while working a little over than a single thread. This implementation is not multi-threaded so the performance is expected to be similar in other invocations. Here again, SNIKT by default, provides an option to speed this step further with looking at top 10,000 reads, in case the user is making an assessment for a significantly larger set. This can be configured by the user as well with the `-s` flag.
- All graphs presented in the report are well labelled and annotated.

<br>

## Comparison of trimming
### NanoFilt execution
The filtering criteria is set to remove any reads below 500 bp and trim by 100 bp from 5' and 20 bp from 3'.
NanoFilt does not read in gzipped fastq files, so a prior decompression is needed.  Since this makes NanoFilt a piped command, it was encapsulated for execution in a script so that the `/usr/bin/time` utility can capture its usage.
```_
$ cat NanoFilt_run.sh
# Add conda init to subshell
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
        . "$(conda info --base)/etc/profile.d/conda.sh"
fi

conda activate nanopack
gunzip -c SRR13441294.fastq.gz | NanoFilt --length 500 --headcrop 100 --tailcrop 20 > SRR13441294_trim_nanofilt.fastq

$ recordStats ./NanoFilt_run.sh
Command stats:
Proc:   Elapsed Time = 7:12.20, Avg CPU = 120%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 82792KB,
Disk:   In = 0, Out =18966872
Exit Status: 0
```
A filtered fastq file was generated from this process. No report was generated.

#### NanoFilt result summary
- The process of trimming and filtering took about 7.2 minutes.
- The resulting data has no statistics generated with it. This means, additional steps might be needed to add a validation or statistical check for the resulting reads. This might add more time and resource usage.

### SNIKT execution
The filtering criteria is set to remove any reads below 500 bp and trim by 100 bp from 5' and 20 bp from 3'.
SNIKT defaults to making a report using pre- and post-trim fastq data. These steps become significant in fraction of time and resources used for the procedure. So, here we ran SNIKT with verbosity that shows time taken by each step.
```_
$ recordStats ./snikt.R -v --filter=500 --trim5=100 --trim3=20 -o SRR13441294_trim_snikt SRR13441294.fastq.gz
[1s] Libraries loaded and user parameters set
[2s] Executing pre-trim seqtk forward pass ... done [2s]
[5s] Executing pre-trim seqtk reverse pass ... done [3s]
[8s] Preparing graphs pre-trimming ... done [1s]
[9s] Executing serial trim and filter command ...
[74s] Trimming finished with status: 0 [65s]
SNIKT cleaned reads are available in:  /home/pr/WD/WD2/SRR13441294_trim_snikt.fastq
[74s] Executing post-trim seqtk forward pass ... done [71s]
[146s] Executing post-trim seqtk reverse pass ... done [67s]
[213s] Preparing graphs post-trimming ... done [50s]
[263s] Rendering report ... rendered [0s]
[263s] Removing intermediate files ... done [0s]
SNIKT contamination report is available in: /home/pr/WD/WD2/SRR13441294_trim_snikt.html

Command stats:
Proc:   Elapsed Time = 4:24.09, Avg CPU = 131%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 764692KB,
Disk:   In = 0, Out =19097272
Exit Status: 0
```
A filtered fastq as well as the following trimming report was generated from this process. Download locally to see.  
[SNIKT post-trim report: ./\_attach/SRR13441294_trim_snikt.html](_attach/SRR13441294_trim_snikt.html)

#### SNIKT result summary
- The entire process, including the report generation steps, finished in about 4.4 minutes. This is significantly faster than 7.2 minutes taken by NanoFilt.
	- If comparing only the trimming and filtering steps, the difference is even larger with SNIKT doing the same task in 65 seconds. This is due to the fact that SNIKT leverages `seqtk` in the background which is an extremely efficient fastq processor.
	- SNIKT still took about 760 MB in RAM memory (NanoFilt took 82 MB) because it is making a post trimming report.
- The trimming report presents a clear summary of the process with demonstration of statistics pre- and post-trimming. No additional steps are needed for validation, simplifying the process.

### Fastq stat comparison post-trimming
Fastq statistics for files generated post-trimming with NanoFilt and SNIKT are summarized below with a custom script.
```_
$ summarizeFastq.pl SRR13441294_trim_*.fastq
#Dataset        #Reads  #Bases  #MinLen #MaxLen #AvgLen #AvgQ   #AvgA   #AvgC   #AvgG   #AvgT   #AvgN
SRR13441294_trim_nanofilt.fastq 398605  4846925580      500     166358  12159.72        21.3    32.8    16.8    17.2    33.1    0.0
SRR13441294_trim_snikt.fastq    398605  4846925580      500     166358  12159.72        21.3    32.8    16.8    17.2    33.1    0.0
```

| Dataset                         | Reads  | Bases      | MinLen | MaxLen | AvgLen   | AvgQ | AvgA | AvgC | AvgG | AvgT | AvgN |
| ------------------------------- | ------ | ---------- | ------ | ------ | -------- | ---- | ---- | ---- | ---- | ---- | ---- |
| SRR13441294_trim_nanofilt.fastq | 398605 | 4846925580 | 500    | 166358 | 12159.72 | 21.3 | 32.8 | 16.8 | 17.2 | 33.1 | 0.0  |
| SRR13441294_trim_snikt.fastq    | 398605 | 4846925580 | 500    | 166358 | 12159.72 | 21.3 | 32.8 | 16.8 | 17.2 | 33.1 | 0.0  |

The results show that both processes - NanoFilt and SNIKT produce fastq reads with the same basic statistics when used the same way, adding a validation to both processes.
