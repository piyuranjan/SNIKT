# Comparison w Nanopack tools
Created: 2022-04-08
Author: Piyush Ranjan

---

Tmp Location: `WD2` (`\\wsl$\Ubuntu\home\pr\WD\\wsl$\Ubuntu\home\pr\WD\WD2`)
Location: `$PWD/Nanopack_Comparison/`

In response to reviewer comments, we compared the performance of SNIKT with two popular tools for Nanopore QC - [NanoQC](https://github.com/wdecoster/nanoQC) and [NanoFilt](https://github.com/wdecoster/nanofilt) that are from the [NanoPack](https://github.com/wdecoster/nanopack) ([DOI](https://doi.org/10.1093/bioinformatics/bty149)) suite of tools.

## Test setup
### Dataset
The dataset used for this comparison is the same SRA run that is used in the main paper. This sequencing run is a WGS long-read Nanopore run with 9.4.1 pore chemistry for HMW DNA extraction from Candida albicans CHN1 which was prepared with the SQK-RAD004 rapid DNA library preparation.
```_
$ summarizeFastq.pl SRR13441294.fastq.gz
#Dataset        #Reads  #Bases  #MinLen #MaxLen #AvgLen #AvgQ   #AvgA   #AvgC   #AvgG   #AvgT   #AvgN
SRR13441294.fastq.gz    412984  4900224854      1       166478  11865.41        21.3    32.7    16.9    17.2    33.2    0.0
```
| Dataset              | Reads  | Bases      | MinLen | MaxLen | AvgLen   | AvgQ | AvgA | AvgC | AvgG | AvgT | AvgN |
| -------------------- | ------ | ---------- | ------ | ------ | -------- | ---- | ---- | ---- | ---- | ---- | ---- |
| SRR13441294.fastq.gz | 412984 | 4900224854 | 1      | 166478 | 11865.41 | 21.3 | 32.7 | 16.9 | 17.2 | 33.2 | 0.0  |

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


## Comparison with NanoQC
### NanoQC Default `--minlen`
Executing NanoQC without `--minlen` as the description seems to suggest it will filter reads when making a QC report.
```_
$ recordStats nanoQC SRR13441294.fastq.gz
BokehDeprecationWarning: plot_width and plot_height was deprecated in Bokeh 2.4.0 and will be removed, use width or height instead.

Command stats:
Proc:   Elapsed Time = 18:31.92,        Avg CPU = 100%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 1759732KB,
Disk:   In = 0, Out =144
Exit Status: 0
```
Following is the report that was generated.
![](_attach/nanoQC.html)

Results suggest the following.
- The default for `--minlen` was taken as 200 and the method did filter reads smaller than that. It made graphs for a 100 bp window from each end. This means that the package necessitates the coupling between the viewing window and the read length filter removing a true representation of the original file before plotting it. The log suggested that it used 411060 reads instead of the original 412984 which leaves 1924 reads outside of the initial representation.
- The graphing step took about 18.5 minutes working over a single thread. This application is not multi threaded (neither is SNIKT).
- While addressable, some graphs are not labelled.
	- While the graphs are interactive, the data drawn on them is fixed. It limits the interactivity to only zooming in. At the same time, zooming in happens at both axis simultaneously, which defeats the purpose of zooming in entirely since the y-axis needs to be scaled to understand the graph.
	- The frequency length graph is set to have the same axis in SNIKT, which allows the user to make a simultaneous decision on how much throughput they would waste at different levels of trimming along with choosing an appropriate criteria using the nucleotide compositions and the Phred scores.


### SNIKT `--notrim`
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
Following was the report that was generated.
![](_attach/SRR13441294_notrim.html)

Results suggest the following.
- Unlike NanoQC, SNIKT took the entire fastq into consideration without any filters. The default x-axis sizes are set for long-read shotgun sequencing, so they work very well in showing the poor quality bases and adapter contamination in an easy view, so that the user can determine where to make their trims. The left and the right panels also demonstrate negligible loss in number of reads (using read frequency) when making a trim at, what we see appropriate, 100 bp from 5' and 20 bp from 3' for this set. The graph is a vector image that can be opened in a bigger view if desired in a new tab or can be saved to zoom in. However, the zoom in will have similar limitations as NanoQC as both axes will zoom together.
- The graphing step took about 4.3 minutes to complete working a little over than a single thread. This implementation is not multi-threaded so the performance is expected to be similar in other invocations. Here again, SNIKT by default, provides an option to speed this step further with looking at top 10,000 reads, in case the user is making an assessment for a significantly larger set. This can be configured by the user as well with the `-s` flag.
- All graphs presented in the report are well labelled and annotated.


## Comparison with NanoFilt
### NanoFilt
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

Results suggest the following.
- The process of trimming and filtering took about 7.2 minutes.
- The resulting data has no statistics generated with is. This means, additional steps might be needed to add a validation or statistical check for the resulting reads. This might add more time and resource usage.

### SNIKT
The filtering criteria is set to remove any reads below 500 bp and trim by 100 bp from 5' and 20 bp from 3'.
SNIKT defaults to making a report using pre- and post-trim fastq data. These steps become significant in fraction of time and resources used for the procedure. So, here we ran SNIKT with verbosity that shows time taken by each step.
```_
$ recordStats ./snikt.R --filter=500 --trim5=100 --trim3=20 -o SRR13441294_trim_snikt SRR13441294.fastq.gz
SNIKT cleaned reads are available in:  /home/pr/WD/WD2/SRR13441294_trim_snikt.fastq
SNIKT contamination report is available in: /home/pr/WD/WD2/SRR13441294_trim_snikt.html

Command stats:
Proc:   Elapsed Time = 4:14.37, Avg CPU = 131%,
Mem:    Avg Total Mem = 0KB,    Peak Mem = 759584KB,
Disk:   In = 119280,    Out =19097272
Exit Status: 0
```

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
A filtered fastq as well as the following trimming report was generated from this process.
![](_attach/SRR13441294_trim_snikt.html)

Results suggest the following.
- The entire process, including the report generation steps, finished in about 4.4 minutes. This is significantly faster than 7.2 minutes taken by NanoFilt.
	- If comparing only the trimming and filtering steps, the difference is even larger with SNIKT doing the same task in 65 seconds. This is due to the fact that SNIKT leverages `seqtk` in the background which is an extremely efficient fastq processor.
- The trimming report presents a clear summary of the process with demonstration of statistics pre- and post-trimming. No additional steps are needed for validation simplifying the process.

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
