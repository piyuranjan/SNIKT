#!/usr/bin/env Rscript

##############
### Config ###
##############

time0 <- Sys.time()

# Library load without warning messages
packages <- c("tidyverse","grid","gridExtra","docopt","lubridate")
# tmp <- lapply(packages,function(x) suppressPackageStartupMessages(require(x,character.only=T)))
tmp <- lapply(packages,function(x) suppressPackageStartupMessages(library(x,character.only=T)))
# tmp <- lapply(packages,function(x) suppressWarnings(suppressMessages(library(x,character.only=T))))
# print(tmp)

# Find exactly how script is executed
initialOptions <- commandArgs(trailingOnly=FALSE)
fileArgName <- "--file="
scriptPath <- sub(fileArgName,"",initialOptions[grep(fileArgName,initialOptions)])
# cat(scriptPath)


#################
### Functions ###
#################

AreaPlot <- function(df, zoomLen=NULL, trim=NULL, pad=0.025, xBreaks=6){
	# Generate an area plot for given compositions per position.
	#   Process either the full length or a zoom length (smaller x-axis) plot.
	#   For full length plots, include average Phred base quality at each position and reads seen
	#     per position (read length distribution) as lines.
	#   Subset compositons by position to process zoom plots and add optional trim lines.
	# Parameters:
	#   df(dataFrame): A data frame with per positon nucleotide counts, nucleotide compositions
	#                  in [ATGCN] or [(AT)(CG)(N)] and Phred scores.
	#   zoomLen(int): Length in NT to subset and fix x-axis positions. 3' or 5' determined later.
	#   trim(int): Position at which to put a trim (vertical dashed) line suggesting trimming location.
	#   pad(float): Left or right plot padding.
	#   xBreaks(int): Number of breaks or ticks on x-axis for the plot.
	# Returns:
	#   ap(ggplot): A ggplot object containing elements for the composition plot.
	
	if(!is.null(zoomLen)){ #subsetting condition to follow for zoom plots
		df <- dplyr::filter(df,POS<=zoomLen)
	}
	
	ap <- ggplot(df,aes(x=POS,y=composition,fill=nt)) +
		geom_area(alpha=0.75) + #area plot; lowered alpha to see gridlines
		# geom_line(aes(y=avgQ),color="black",alpha=0.75) + #avg Phred score
		geom_line(aes(y=errQ),color="black",alpha=0.75) +  # use errQ instead of avgQ as better Q estimator
		geom_line(aes(y=bases/max(bases)*100),color="red") + #make normalized length freq as line
		scale_y_continuous(labels=c(0,25,50,75,100),expand=c(pad,0), #adjust top/bottom plot padding
			sec.axis=sec_axis(~.*max(df$bases)/100,name="Frequency",labels=FormatSI())) + #secondary axis uses SI scale
		# theme_bw() + #remove gray background
		theme_linedraw() +
		theme(legend.position="none",title=element_blank(),axis.line=element_line(color="black"), #remove legend & titles; set axis colors
			axis.line.y.right=element_line(color="red"),axis.ticks.y.right=element_line(color="red"),
			axis.text.y.right=element_text(color="red"))
	
	# Set fill colors and determine x-axis length
	if(length(unique(df$nt))==5){  # condition for 4NT+1
		ap <- ap + scale_fill_manual(
			values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")  # fill order: "A", "G", "C", "T", "N"
		)
		maxX <- dim(df)[1]/5  # length of x-axis to which data is present for 4NT+1
	} else{  # assume 2NT+1 (length == 3)
		ap <- ap + scale_fill_manual(values=c("#00274C","#FFCB05","#E76BF3"))  # fill order: "AT", "GC", "N"
		maxX <- dim(df)[1]/3  # length of x-axis to which data is present for 2NT+1
	}
	
	# Condition to skip for zoom plots, return center plot.
	if(is.null(zoomLen)){
		ap <- ap + scale_x_continuous(
			expand=c(pad,0),  # adjust left/right plot padding
			breaks=BreaksWMax(scales::breaks_pretty(n=xBreaks),maxX)  # adjust x-axis breaks
		)
		return(ap) #return center plot
	}
	
	# Positioned here to avoid null zoomLen comparison.
	# From zoom argument, determine if forward or reverse scale applies
	if(deparse(substitute(zoomLen))=="zoom5Len")
		{ap <- ap + scale_x_continuous(expand=c(pad,0), breaks=scales::breaks_pretty(n=xBreaks))} #adjust left/right plot padding
	else #assumes zoom3Len
		{ap <- ap + scale_x_reverse(expand=c(pad,0), breaks=scales::breaks_pretty(n=xBreaks))} #reverse x-axis, adjust left/right plot padding
	# Add dashed vertical line at trim position
	if(!is.null(trim)){ap <- ap + geom_vline(xintercept=trim,linetype="dashed")}
	return(ap) #return zoom plot
}

BreaksWMax <- function(originalFun = scales::breaks_pretty(), maxX){
	# Add the max value of a graph to breaks.
	# Parameters:
	#   originalFun(function): a function call to scales::breaks_pretty()
	#   maxX(int): max break on x-axis to put a tick
	# Returns:
	#   breaksSort: breaks in sorted order with a max value
	
	function(x){
		originalBreaks <- originalFun(x)
		# breaks <- c(originalBreaks, as.integer(max(x)))  # this max does not correspond to data on x-axis
		breaks <- c(originalBreaks, as.integer(maxX))  # add max value to existing breaks
		breaksSort <- sort(breaks)
		# Remove any breaks that are less than 50 bp apart
		close <- diff(breaksSort) < 50
		breaksSort <- breaksSort[!close]
		return(breaksSort)
	}
}

CreateGrid <- function(
	legendObj=legends, p1L=areaComp4Zoom5, p1R=areaComp4Zoom3, p2L=areaComp2Zoom5,
	p2R=areaComp2Zoom3, p1C=NULL, p2C=NULL, centerLen=NA){
	# Create and return a grid of 4 or 6 plots with axis labels and legends.
	#   This step is the most time consuming as plot objects are finalized to pixels.
	#   No change in plot objects is possible after execution of this function.
	# Parameters:
	#   legendObj(ggplotLegend): Legend part stripped from a ggplot object to serve as
	#                            legends for the figure constructed here.
	#   p1L(ggplot): A ggplot object containing 5' zoom plot positions for 4 NTs(A,T,G,C) + N.
	#   p1R(ggplot): A ggplot object containing 3' zoom plot positions for 4 NTs(A,T,G,C) + N.
	#   p2L(ggplot): A ggplot object containing 5' zoom plot positions for 2 NT combinations(AT,GC) + N.
	#   p2R(ggplot): A ggplot object containing 3' zoom plot positions for 2 NT combinations(AT,GC) + N.
	#   p1C(ggplot): A ggplot object containing all positions for 4 NTs(A,T,G,C) + N.
	#   p2C(ggplot): A ggplot object containing all positions for 2 NT combinations(AT,GC) + N.
	#   centerLen(int): Length of the x-axis in center plot that has all positions.
	# Returns:
	#   (grid.arrange): Aligned grid with 4 or 6 plots, attached axis labels and legends.
	
	# Make labels for the plot grid
	topLabel <- textGrob(label="Aggregate Sequence Representation Per Nucleotide Position",
		gp=gpar(fontface=2,fontsize=20))
	bottomLabel <- textGrob(label="Read length (nucleotide positions)",
		gp=gpar(fontface=2,fontsize=14))
	leftLabel <- textGrob(label="Nucleotide composition (0-100), Average Phred score",
		rot=90,gp=gpar(fontface=2,fontsize=14))
	rightLabel <- textGrob(label="Read frequency (reads available at each position)",
		rot=90,gp=gpar(fontface=2,fontsize=14,col="red"))
	
	# Make left and right plot groups
	leftPlots <- arrangeGrob(p1L,p2L,ncol=1,
		top=textGrob("Aligned 5' beginning",gp=gpar(fontface=2,fontsize=14)))
	rightPlots <- arrangeGrob(p1R,p2R,ncol=1,
		top=textGrob("Aligned 3' ending",gp=gpar(fontface=2,fontsize=14)))
	
	if(is.null(p1C)){ #bind 4 plots in grid
		areaArranged <- arrangeGrob(leftPlots,rightPlots,ncol=2,
			top=topLabel,bottom=bottomLabel,left=leftLabel,right=rightLabel)
	} else{
		# Make center plot groups and bind 6 plots in grid
		centerPlots <- arrangeGrob(p1C,p2C,ncol=1,
			top=textGrob(paste("Aligned 5' to length",centerLen),gp=gpar(fontface=2,fontsize=14)))
		areaArranged <- arrangeGrob(leftPlots,centerPlots,rightPlots,ncol=3,widths=c(1,2,1),
			top=topLabel,bottom=bottomLabel,left=leftLabel,right=rightLabel)
	}
	# Render 4 or 6 plot grid and return
	# This step takes the most time
	grid.arrange(areaArranged,legendObj,nrow=2,heights=c(28,2))
}

CreateLegend <- function(){
	# Create a legend object that can be passed to all plots.
	# Parameters:
	#   None
	# Returns:
	#   (object): A ggplot legend object
	
	# Make a dummy plot with display settings
	g <- ggplot(
			data.frame(NT=fct_inorder(paste0("% ",c("A","G","C","T","AT","GC","N"))),N=1:7),
			aes(x=NT,y=N,fill=NT)
		) + 
		geom_col(alpha=0.75) + 
		geom_line(aes(group=1,color="Average Phred Score  ")) +
		geom_line(aes(y=(N-.5),group=1,color="Nucleotide Per Position")) + 
		theme_bw() +
		scale_color_manual(
			name=element_blank(),
			breaks=c("Average Phred Score  ","Nucleotide Per Position"),
			values=c("Average Phred Score  "="black","Nucleotide Per Position"="red")
		) + #create line legend
		scale_fill_manual(
			values=c("#F8766D","#A3A500","#00BF7D","#00B0F6","#00274C","#FFCB05","#E76BF3")  # "A", "G", "C", "T", "AT", "GC", "N" 
		) +
		guides(fill=guide_legend(title=NULL,order=1,nrow=1),color=guide_legend(order=2)) +
		theme(legend.position="bottom")
	# Extract the legend and return
	leg_grob <- ggplotGrob(g)$grobs
	leg_grob[[which(map_chr(leg_grob,function(x) x$name)=="guide-box")]]
}

#--- Needs Review ---#
# Likely not needed anymore
DecompressGzip <- function(fileName, path=tmpDir, prefix=outPrefix, proc=1){
	# Decompress a gzipped fastq file and return new fileName
	# Parameters:
	#   fileName(str): Gzipped fastq file (w/ path) to process
	#   path(str): Path where decompressed fastq should be kept
	#   prefix(str): Prefix for naming the decompressed fastq
	#   proc(int): Number of threads to use for decompression
	# Returns:
	#   outFile(str): Decompressed fastq fileName w/ path
	
	outFile <- paste0(path,"/",prefix,".fastq")
	pigzEC <- system("pigz -V",ignore.stderr=TRUE) #check if pigz is available
	if(pigzEC==0){
		cmd <- paste("pigz -p",proc,"-dc",fileName,">",outFile) #use pigz or
	} else{
		cmd <- paste("gzip -dc",fileName,">",outFile) #use gzip
	}
	# cat(paste0("cmd: ",cmd,"\n"))
	cmdEC <- system(cmd)
	if(cmdEC!=0){stop(paste0("ERR: Decompression exited with code: ",cmdEC))}
	return(outFile)
}
#--- ---#

ExecSeqtk <- function(fastq=seqFile, revComp=FALSE, head=0, isFastqGz=isFqGz){
	# Run seqtk on a file and return buffered output.
	# Parameters:
	#   fastq(str): Fastq filename, with path, to process on.
	#   revComp(logi): Should the sequences be reverse complemented before feeding to seqtk?
	#   head(int): Number of fastq records to pass to seqtk.
	#   isFastqGz(logi): Is the incoming fastq gzipped?
	# Returns:
	#   seqtkBuffer(str): Output from seqtk as a buffered string.
	
	cmd=""
	fastq <- str_replace_all(fastq,' ','\\\\ ') # escape any spaces in file or path
	# writeLines(paste('fastq:',fastq))
	# Prepare command to handle taking top fastq records and decompression
	if(head>0){ #decompression must be done with head
		lines <- head*4 #convert fastq records to lines in file
		if(isFastqGz){
			cmd <- paste0(cmd,"gzip -dc ",fastq," | head -",lines," | ")
		} else{
			cmd <- paste0(cmd,"head -",lines," ",fastq," | ")
		}
	} else{ #decompression will be done by seqtk if full fastq is being processed
		cmd <- paste0(cmd,"cat ",fastq," | ")
	}
	# writeLines(paste('cmd:',cmd))
	# q('no',0,FALSE)
	# Prepare command for seqtk
	if(revComp==FALSE){ #compute forward 5' aligned compositions
		cmd <- paste0(cmd,"seqtk fqchk -")
	} else{ #compute reverse complement for 3' alignment and calculate compositions
		cmd <- paste0(cmd,"seqtk seq -r - | seqtk fqchk -")
	}
	# cat(cmd)
	seqtkBuffer <- system(cmd,intern=T) #execute seqtk on shell and retrieve data
	if(debug){write(paste0('\nDEBUG: ',cmd),stderr())}
	return(seqtkBuffer)
}

ExtractSummary <- function(buffer=seqtkBuffer){
	# Extract summary statistics from a seqtk buffer and return it.
	# Parameters:
	#   buffer(str): Buffered output from a seqtk fqchk run on a fastq file
	# Returns:
	#   summDf(dataFrame): A data frame with all summary values organized as needed
	
	# Check seqtk output format consistency
	colOrder <- "POS	#bases	%A	%C	%G	%T	%N	avgQ	errQ	%low	%high"
	if(buffer[2]!=colOrder){stop("ERR: seqtk output looks inconsistent. Please check seqtk version.")}
	# Extract specific fields from initial lines of the seqtk buffer
	if(length(str_extract_all(buffer[1],"(?<=:\\s)\\d+\\.?\\d*")[[1]]) == 3){
		summ1 <- str_extract_all(buffer[1],"(?<=:\\s)\\d+\\.?\\d*")[[1]] %>%
			setNames(c("Min Length","Max Length","Avg Length"))
		# summ2 <- str_split(buffer[3],"\t")[[1]][2:8] %>%
		summ2 <- str_split(buffer[3],"\t")[[1]][c(2,3,4,5,6,7,9)] %>%  # read errQ as average quality
			setNames(c("Num Bases","A","C","G","T","N","Avg Q"))
		summ3 <- str_split(buffer[4],"\t")[[1]][2] %>%
			setNames("Num Seq")
	} else{
		# write(paste0("ERR: This file is either not a fastq or doesn't have any reads.\n",stderr()))
		stop("ERR: This file is either not a fastq or doesn't have any reads.")
	}
	# summ1 <- str_extract_all(buffer[1],"(?<=:\\s)\\d+\\.?\\d*")[[1]] %>%
	# 	setNames(c("Min Length","Max Length","Avg Length"))
	# summ2 <- str_split(buffer[3],"\t")[[1]][2:8] %>%
	# 	setNames(c("Num Bases","A","C","G","T","N","Avg Q"))
	# summ3 <- str_split(buffer[4],"\t")[[1]][2] %>%
	# 	setNames("Num Seq")
	# Collect all summary values and order them as needed
	summDf <- t(c(summ3,summ1,summ2)) %>%
		data.frame(check.names=F,stringsAsFactors=F) %>%
		select(1,5,2:4,11,everything())
	# print(summDf)
	return(summDf)
}

ExtractCompositions <- function(buffer=seqtkBuffer, fractionHide=percHide){
	# Extract nucleotide compositions from a seqtk buffer and return it.
	#   Also features removing super long, low occupancy 3' positions making
	#   graph rendering significantly faster.
	# Parameters:
	#   buffer(str): Buffered output from a seqtk fqchk run on a fastq file
	#   fractionHide(float): Fraction of total bases by which to filter the long end tail
	# Returns:
	#   df(dataFrame): A data frame with positional compositions for 4 NT and 2 NT combinations
	
	# Process the seqtk buffer, store and subset
	colNames <- c("POS","bases","A","C","G","T","N","avgQ","errQ","low","high") #full list of columns in seqtk output
	df <- read_delim(I(buffer),"\t",skip=3,col_names=colNames,show_col_types = FALSE) #read leaving first 3 summary lines
	# df <- df[,1:8] #leave out last 3 unneccesary columns
	df <- df[,1:9] #leave out last 2 unneccesary columns
	
	# Filter compositions by percentage of total bases
	subdf <- df[nrow(df):1,1:2] %>% #subset cols 1,2 and reverse rows
		mutate(CumBases=cumsum(bases)) #calcualte cumulative sum of bases
	numFilter <- (subdf[nrow(subdf),3]*fractionHide)[[1]] #calculate bases under the curve for filter threshold
	trimPoint <- nrow(subdf[subdf[,3]>numFilter,])
	df <- dplyr::filter(df,POS<=trimPoint) #remove low confidence compositions by filtering the long tail
	
	df <- mutate(df,AT=A+T,GC=G+C) #add paired AT and GC compositions
	return(df)
}

FormatSI <- function(...) {
	# source: https://www.moeding.net/2011/10/metric-prefixes-for-ggplot2-scales/
	# Format a vector of numeric values according to the International System of Units.
	# http://en.wikipedia.org/wiki/SI_prefix
	#
	# Based on code by Ben Tupper: https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
	# Args: ...: Args passed to format()
	#
	# Returns: A function to format a vector of strings using SI prefix notation
	
	function(x) {
		limits <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e0, 1e3, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21, 1e24)
		prefix <- c("y",   "z",   "a",   "f",   "p",   "n",  "µ",  "m",  "",  "K", "M", "G", "T",  "P",  "E",  "Z",  "Y")
		
		# Vector with array indices according to position in intervals
		i <- findInterval(abs(x),limits)
		
		# Set prefix to " " for very small values < 1e-24
		i <- ifelse(i==0,which(limits==1e0),i)
		paste(format(round(x/limits[i],1),trim=TRUE,scientific=FALSE,...),prefix[i])
	}
}

ReadTrimLen <- function(msg="Enter the trimming length [>=0]: ", maxLen){
	# Read a whole number in a range from user input and return it.
	# Parameters:
	#   msg(str): Prompt message for the user. Embedding possible range in this msg is recommended.
	#   maxLen(int): Max possible range for the user input.
	# Returns:
	#   n(int): A positive whole number.
	
	cat(msg)
	con <- file("stdin")
	n <- readLines(con, n=1)
	if(!grepl("^[0-9]+$",n)){ #check for +ve integers
		write("WARN: Enter a whole number [>=0] please.\n",stderr())
		n <- ReadTrimLen(msg,maxLen)
	} else if(as.integer(n)>maxLen){ #check for trim to be shorter than max possible
		write(paste0("WARN: Trim length can not be longer than: ",maxLen,"\n"),stderr())
		cat("One of the two scenarios is happening with the provided trim length.\n")
		cat("1 - It's larger so it overlaps the trim on the other side.\n")
		cat("2 - It's longer than the zoom length which suggests that you may not be using the temporary graph for inference.\n")
		cat("If you want to use your own trim length, please run the program in non-interactive mode by specifying -t5, -t3 at execution. SNIKT will adjust zoom lengths accordingly.\n")
		n <- ReadTrimLen(msg,maxLen)
	}
	close(con)
	return(as.integer(n))
}

RenderMd <- function(conn=mdConn, mdtxt=mdFile, mdhtml=mdReport, tmp=tmpDir, keepTmp=keepTmpDir){
	# Render the markdown report, clear intermediates, finish with a success code.
	# Parameters:
	#   conn(fileHandle): File handle that writes to the text markdown file
	#   mdtxt(str): File name w/ path for the text markdown
	#   mdhtml(str): File name w/ path for the HTML report
	#   tmp(str): Path for the temporary directory storing any intermediates
	#   keepTmp(logi): Should temporary directory be kept?
	# Returns:
	#   None: Quits the program with a successful execution code 0.
	
	writeLines(paste0("---\n\n<br>\n\nTo see figures in larger view, open in new tab or save."),conn)
	writeLines(paste0("\n\nReport created: ",Sys.time()),conn)
	writeLines(paste0("\n\n",programVersion),conn)
	close(conn) #close the connection to text markdown
	if(verbose) cat(TimeDiff(time0),"Rendering report ... ")
	timeStep<-Sys.time()
	html <- rmarkdown::render(mdtxt,output_file=mdhtml,quiet=TRUE) #render the text markdown to HTML
	if(verbose) cat("rendered",TimeDiff(timeStep),"\n")
	
	# Remove extra Rplots.pdf if generated by this instance; relates to a bug
	# No need for this resolution if tmpDir is set as current working dir
	# if(file.exists("Rplots.pdf")){invisible(file.remove("Rplots.pdf"))}
		
	# Remove temporary dir and files
	if(!keepTmp){
		if(verbose) cat(TimeDiff(time0),"Removing intermediate files ... ")
		timeStep <- Sys.time()
		setwd(dirname(mdhtml)) #jump cwd to HTML report path before deleting tmp
		unlinkStat <- unlink(tmp,recursive=TRUE)
		if(unlinkStat==0){
			if(verbose) cat("done",TimeDiff(timeStep),"\n")
		} else{
			write(paste0("\nWARN: Could not delete intermediates: ",tmp),stderr())
		}
	} else{
		if(verbose) cat("Intermediate files are available in:",tmp,"\n")
	}
	
	cat("SNIKT contamination report is available in:",html,"\n")
	q("no",0,FALSE) #quit the program saying successful execution
}

TimeDiff <- function(t0, t1=Sys.time()){
	# Return a formatted time difference between two time points.
	# Parameters:
	#   t0(Sys.time): A Sys.time() object serving as the starting point.
	#   t1(Sys.time): A Sys.time() object serving as the ending point.
	# Returns:
	#   diffStamp(str): A formatted string with time in seconds.
	
	if(is.null(t0)){ #quit with error if t0 is not provided
		write(paste0("ERR: Need t0 to calculate time diff.\n"),stderr())
		q("no",1,FALSE)
	}
	diffStamp <- paste0("[",round(time_length(t1-t0,unit="second"),0),"s]")
	return(diffStamp)
}

TimeStamp <- function(currTime=Sys.time()){
	# Return a formatted time stamp w/o any special chars.
	# Parameters:
	#   currTime(Sys.time): A Sys.time() object to convert, default to current time.
	# Returns:
	#   (str): A formatted time stamp in format yyyymmddhhmmss.
	
	gsub("[: -]", "", currTime, perl=TRUE)
}

TrimFiltSeq <- function(seqFile, trim5Len, trim3Len, filtLen, trimFile=NULL){
	# Trim and filter fastq reads with seqtk trimfq and seqtk seq.
	# Parameters:
	#   seqFile(str): Fastq file name w/ path to trim.
	#   trim5Len(int): Nucleotide length to trim from 5' end.
	#   trim3Len(int): Nucleotide length to trim from 3' end.
	#   filtLen(int): Reads less than this length will be dropped.
	#   trimFile(str): Optional file name w/path for the trimmed file.
	# Returns:
	#   trimFile(str): Trimmed fastq file name w/ path.
	
	# Creare trimmed file name if not passed during function call
	if(is.null(trimFile)){trimFile <- sub("\\.f[ast]{0,3}(a|q)(\\.gz)?","_trim.fastq",seqFile,perl=TRUE)}
	
	# Make trimming and filtering command
	fastq <- str_replace_all(trimFile,' ','\\\\ ') # escape any spaces in file or path
	seqFile <- str_replace_all(seqFile,' ','\\\\ ') # escape any spaces in file or path
	cmd <- paste0("seqtk trimfq -b ",trim5Len," -e ",trim3Len," ",seqFile)
	statusMsg <- paste0("Executing serial trim ")
	if(filtLen>0){
		cmd <- paste0(cmd," | seqtk seq -L ",filtLen," -")
		statusMsg <- paste0(statusMsg,"and filter ")
	}
	cmd <- paste0(cmd," >",fastq)
	statusMsg <- paste0(statusMsg,"command ...")
	if(debug){statusMsg <- paste0(statusMsg,"\nDEBUG: ",cmd)} #embed seqtk command in status message
	
	if(verbose) cat(TimeDiff(time0),statusMsg,"\n")
	timeStep<-Sys.time()
	
	# Execute and quit with error if seqtk was unsuccessful
	exitCode <- system(cmd)
	if(verbose) cat(TimeDiff(time0),"Trimming finished with status:",exitCode,TimeDiff(timeStep),"\n")
	if(exitCode>0){ 
		write(paste0("ERR: seqtk trimming unsucessful; seqtk exit code: ",exitCode,"\n"),stderr())
		q("no",1,FALSE)
	}
	
	return(trimFile)
}

#--- Needs review ---#
# Chaining trimming commands with GNU Parallel was significantly less efficient than reads
#   directly being handled by seqtk in a single thread.
# There is no need for this function and so is deprecated.
TrimSeqParallel <- function(seqFile, trim5Len, trim3Len, proc=2, fqPerThread=100){
	# Trim and filter fastq reads with seqtk trimfq and seqtk seq on parallel threads.
	#   Uses GNU Parallel for thread management in Unix.
	# Parameters:
	#   seqFile(str): Fastq file name w/ path to trim.
	#   trim5Len(int): Nucleotide length to trim from 5' end.
	#   trim3Len(int): Nucleotide length to trim from 3' end.
	#   proc(int): Number of threads to use.
	#   fqPerThread(int): Number of fastq records to pass on one thread.
	# Returns:
	#   trimFile(str): Trimmed fastq file name w/ path.
	
	# Prepare command
	trimFile <- sub("\\.f[ast]{0,3}(a|q)(\\.gz)?","-trim.fastq",seqFile,perl=T)
	cmd <- ""
	
	# Handle decompression
	if(grepl("\\.gz$",seqFile,perl=T)){
		cmd <- paste0(cmd,"gzip -dc ",seqFile)
	} else{
		cmd <- paste0(cmd,"cat ",seqFile)
	}
	
	seqtkCmd <- paste0("seqtk trimfq -b ",trim5Len," -e ",trim3Len," -") #seqtk command with stdin
	
	# Handle parts of GNU Parallel
	fqLines <- fqPerThread*4
	cmd <- paste0(cmd,"|parallel -j ",proc," -k --pipe -N ",fqLines)
	cmd <- paste0(cmd," '",seqtkCmd,"' >",trimFile)
	
	# cat(paste0(cmd,"\n"))
	if(verbose) cat(TimeDiff(time0),"Executing parallel trimming command:",cmd,"\n")
	timeStep<-Sys.time()
	
	# Execute and quit with error if seqtk was unsuccessful
	exitCode <- system(cmd)
	if(verbose) cat(TimeDiff(time0),"Trimming finished with status:",exitCode,TimeDiff(timeStep),"\n")
	if(exitCode>0){
		write(paste0("ERR: seqtk trimming unsucessful; seqtk exit code: ",exitCode,"\n"),stderr())
		q("no",1,FALSE)
	}
	
	return(trimFile)
}
#--- ---#



############
### Main ###
############

# Status message that records time for libraries to load. Feel free to uncomment
#   the following line if you want to see it.
# cat(TimeDiff(time0),"Libraries loaded\n")


## Defaults and user defined parameters ##

## Set docopt string help and read user arguments
docString <- 
"
SNIKT: FastQ QC and sequence over-representation check.
       A wrapper around seqtk to plot per-position nucleotide composition
       for finding and trimming adapter contamination in fastq reads.
       Also filters reads by a length threshold.
Authors Piyush Ranjan, Christopher Brown

For first-time users, interactive mode is recommended.
For detailed help and examples, please visit
https://github.com/piyuranjan/SNIKT

Location: scriptPath

Usage:
  snikt.R [options] [--] <fastq>
  snikt.R <fastq>  # Interactive # Easiest
  snikt.R [--zoom5=<nt>] [--zoom3=<nt>] <fastq>  # Interactive
  snikt.R [(--trim5=<nt> --trim3=<nt>) | --notrim] <fastq>
  snikt.R [--illumina] [-n] <fastq>

Input:
  <fastq>               Sequence file in fastQ format with exts: .fq, .fq.gz,
                          .fastq, .fastq.gz

Options:
  Presets:
  --illumina            This presets options that are better for short-read
                          Illumina datasets.
                          Sets: -f 0 -Z 50 -z 50 --hide=0
                          Defaults are configured for long-read Nanopore fastq.

  Graphing:
  --hide=<frac>         Hide the composition tail by a fraction of total bases.
                          Significantly improves speed, removes end-tail (3')
                          distortion for variable length read sets.
                          [range: 0..1] [default:0.01]
  -s, --skim=<num>      Use top num reads for pre- or no-trim graphs. This
                          improves speed. No effect on post-trim graphs.
                          Use 0 to disable skimming and utilize all reads.
                          [range: 0..maxFastqReads] [default: 10000]
  -x, --xbreaks=<num>   Suggest number of ticks/breaks on x-axis in all graphs.
                          Can be set if the breaks or gridlines are too sparse
                          or dense to determine appropriate trimming. Internal
                          algorithm adjusts ticks.
                          [range: 1..10] [default: 6]
  -Z, --zoom5=<nt>      Zoom-in from aligned 5' beginning to nt bases.
                          [range: 1..maxSeqLen] [default:300]
  -z, --zoom3=<nt>      Zoom-in from aligned 3' ending to nt bases.
                          [range: 1..maxSeqLen] [default:100]
  QC:
  -f, --filter=<nt>     Filter (drop) reads with length < nt after any trimming.
                          [range: 0..maxSeqLen] [default:500]
  -n, --notrim          Disable positional trimming; useful for short-read data
                          Takes precedence over and sets: -T 0 -t 0
  -T, --trim5=<nt>      Trim nt bases from aligned 5' side.
                          [range: 0..(maxSeqLen-trim3)] [default: interactive]
  -t, --trim3=<nt>      Trim nt bases from aligned 3' side.
                          [range: 0..(maxSeqLen-trim5)] [default: interactive]
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
"
# Dynamic modification of docString below Usage section results in usage errors
docString <- sub("scriptPath",scriptPath,docString) #add script location with how it is fired
# cat(docString)
version <- '0.5.0'
programVersion <- paste0("SNIKT ",version,"\n")
arg <- docopt(docString,version=programVersion)
# print(arg)

# Read only verbose and debug from user first. The earlier these are set the better.
verbose <- arg$verbose
debug <- arg$debug
if(debug){verbose <- TRUE}
# Allow traceback for errors if debug is TRUE
if(debug){
	options(error=function(){
		traceback(3)
		if(!interactive()) {q("no",1,FALSE)}
	})
}
# Print docopt parsed values
if(debug){
	write("DEBUG: Docopt values as they came from user input.",stdout())
	print(arg)
}
# q("no",0,FALSE)


## Set defaults before argument parsing
zoom5Len <- as.integer(NA)
zoom3Len <- as.integer(NA)
filterLen <- as.integer(NA)
percHide <- as.numeric(NA)


## Parse arguments in variables

# Required variables
seqFile <- arg$fastq
seqFileName <- basename(seqFile)
# writeLines(paste("seqFile",seqFile,"\nseqFileName",seqFileName))
# q("no",0,FALSE)

# Preset variables, their arguments and DEFAULTS
illumina <- arg$illumina
if(illumina){ #set variables for Illumina reads
	filterLen <- 0
	zoom5Len <- 50
	zoom3Len <- 50
	percHide <- 0
}
if(debug) write(paste0("DEBUG: Values after parsing --illumina:\t-Z=",zoom5Len,"\t-z=",zoom3Len,"\t-f=",filterLen,"\t--hide=",percHide,"\n"),stderr())
# Arguments used in presets need special parsing
if(!is.null(arg$zoom5)){
	zoom5Len <- as.integer(arg$zoom5)
} else if(is.na(zoom5Len)){
	zoom5Len <- 300 #default for zoom5Len (-Z) must be set here. Setting only in docopt string won't work.
}
if(!is.null(arg$zoom3)){
	zoom3Len <- as.integer(arg$zoom3)
} else if(is.na(zoom3Len)){
	zoom3Len <- 100 #default for zoom3Len (-z) must be set here. Setting only in docopt string won't work.
}
if(!is.null(arg$filter)){
	filterLen <- as.integer(arg$filter)
} else if(is.na(filterLen)){
	filterLen <- 500 #default for filterLen (-f) must be set here. Setting only in docopt string won't work.
}
if(!is.null(arg$hide)){ #percent bases to hide from 3'. Avoids rendering noise at the end for full graph.
	percHide <- as.numeric(arg$hide)
} else if(is.na(percHide)){
	percHide <- 0.01 #default for percHide (--hide) must be set here. Setting only in docopt string won't work.
}
if(debug) write(paste0("DEBUG: After parsing preset values:\t-Z=",zoom5Len,"\t-z=",zoom3Len,"\t-f=",filterLen,"\t--hide=",percHide,"\n"),stderr())

# Graphing variables
headFastq <- as.integer(arg$skim)  # use only top (head) reads for initial inference
xBreaks <- as.integer(arg$xbreaks)  # add x-axis breaks/ticks for all graphs

# print(arg)
# q("no",0,FALSE)

# QC variables
noTrim <- arg$notrim
if(arg$notrim){
	trim5Len <- 0; trim3Len <- 0
} else if((arg$trim5=="interactive")|(arg$trim3=="interactive")){
	trim5Len <- as.integer(NA); trim3Len <- as.integer(NA)
} else{
	trim5Len <- as.integer(arg$trim5); trim3Len <- as.integer(arg$trim3)
	if(arg$trim5<0 | arg$trim3<0){stop("ERR: -T or -t out of range.")}
}

# IO variables
if(arg$out=="fastqNoExtension"){
	outPrefix <- sub("\\.f[ast]{0,3}q(\\.gz)?","",seqFileName,perl=TRUE) #extract sequence file prefix
	trimFileName <- paste0(outPrefix,"_trim.fastq") #make trim filename
} else{
	outPrefix <- arg$out
	trimFileName <- paste0(outPrefix,".fastq") #make trim filename
}
if(arg$workdir=="./"){
	workDir <- getwd()
} else{
	workDir <- sub("\\/$","",arg$workdir) #extract path without trailing slash
}
isFqGz <- arg$gzip
tmpDir <- paste0(workDir,"/",outPrefix,"_SniktTemp_",TimeStamp(time0))
keepTmpDir <- arg$keep
# writeLines(paste('workDir',workDir,'\ntmpDir',tmpDir))
# q('no',0,FALSE)

# Generic variables
# threads <- as.integer(arg$proc)


## Defaults and conditions after argument parsing
seqFile <- normalizePath(seqFile) #resolve absolute path
# seqFile <- str_replace_all(seqFile,' ','\\\\ ') # escape any spaces in file or path
if(!file.exists(seqFile)){stop("ERR: Need a sequence file to proceed. See help with -h.")} #validation check for seqFile
if(grepl("\\.gz(ip)?$",seqFile)) isFqGz <- TRUE #override isGzip if extension found in file name
if(percHide>1|percHide<0){stop("ERR: --hide out of range.")}
if(percHide>0.5){write(paste0("WARN: Hiding too much sequence data is not recommended. --hide is currently: ",percHide),stderr())}
if(headFastq<0){stop("ERR: -s out of range, needs to be a positive number.")}
if(xBreaks<0){stop("ERR: -x out of range, needs to be a positive number.")}
if(zoom5Len<1|zoom3Len<1){stop("ERR: -Z or -z out of range.")}
if(((arg$trim5=="interactive")&(arg$trim3!="interactive"))|((arg$trim5!="interactive")&(arg$trim3=="interactive"))){ #warning on only setting one trim length
	write(paste0("WARN: Only one trim length is set. Forcing interactive mode."),stderr())
}
if(filterLen<0){stop("ERR: -f out of range.")}
# if(threads<1){stop("ERR: -p out of range.")}

# Resolve bug regarding generation of Rplots.pdf
# No need for this block if current dir is being set later
# if(file.exists("Rplots.pdf")){ #renaming any preexisting Rplots.pdf so that it doesn't overwrite
	# invisible(file.rename("Rplots.pdf",paste0("Rplots_",TimeStamp(),".pdf")))
# }

if(verbose) cat(TimeDiff(time0),"Libraries loaded and user parameters set\n")
# q('no',0,FALSE)  # user parameters have been set at this point


# Setup workspace variables
dir.create(tmpDir)
# Set work dir to the tmpDir, which will move generation of Rplots.pdf to this folder
# After this, the relative path of every file will change. Use absolute paths for files.
setwd(tmpDir)
legends <- CreateLegend() #extract legends from a simple plot to use in every plot later

#--- Needs review ---#
# # Decompress (extract) if sequence file is compressed
# if(isFqGz){ #if file is known to be gzipped
# 	if(verbose) cat(TimeDiff(time0),"Decompressing fastq ... ")
# 	timeStep <- Sys.time()
# 	seqFile <- DecompressGzip(seqFile,tmpDir,outPrefix,threads)
# 	if(verbose) cat("done",TimeDiff(timeStep),"\n")
# }
#--- ---#


## Execute 1st round of seqtk and gather summary statistics ##

if(verbose) cat(TimeDiff(time0),"Executing pre-trim seqtk forward pass ... ")
timeStep <- Sys.time()
seqtkBuffer <- ExecSeqtk(seqFile,head=headFastq) #compute compositions via seqtk
if(verbose) cat("done",TimeDiff(timeStep),"\n")
# summRawFq <- ExtractSummary(seqtkBuffer,seqFileName) #extract fastq summary from seqtkBuffer
summRawFq <- ExtractSummary(seqtkBuffer) #extract fastq summary from seqtkBuffer
maxSeqLen <- as.integer(summRawFq$`Max Length`)

# Correct zoom and trim size according to summary parameters
# Non-inteactive runs will have trim lengths available at this time.
# Sum of trim legnths cannot exceed max sequence length 
if((!is.na(trim5Len))&(!is.na(trim3Len))&((trim5Len+trim3Len)>=maxSeqLen)){
	stop(paste0(
		"ERR: trim5 (",trim5Len,") length + trim3 (",trim3Len,") length >= max seq length (",maxSeqLen,").\n",
		"Please select non-overlapping trim lengths with -T, -t."
	))
}
# Individual zoom window cannot be larger than the max sequence length
if(zoom5Len>maxSeqLen){
	zoom5Len <- maxSeqLen-1
	write(paste0("WARN: zoom5 larger than the max seq length. Setting -Z to ",zoom5Len,"."),stderr())
}
if(zoom3Len>maxSeqLen){
	zoom3Len <- maxSeqLen-1
	write(paste0("WARN: zoom3 larger than the max seq length. Setting -z to ",zoom3Len,"."),stderr())
}
# Zoom lengths (set or unset by the user) cannot be smaller than the trim lengths
#   in the non-interactive run, for side graphs to appear properly.
if((!is.na(trim5Len))&(zoom5Len<trim5Len)){
	zoom5Len <- trim5Len
	write(paste0("WARN: zoom5 smaller than the specified trim5 length. Setting -Z to ",zoom5Len,"."),stderr())
}
if((!is.na(trim3Len))&(zoom3Len<trim3Len)){
	zoom3Len <- trim3Len
	write(paste0("WARN: zoom3 smaller than the specified trim3 length. Setting -z to ",zoom3Len,"."),stderr())
}

# Compute nucleotide compositions
rawComp <- ExtractCompositions(seqtkBuffer,percHide)  #extract fastq compositions from seqtkBuffer
# rm(seqtkBuffer) #destroy seqtk buffer for memory

# Transform NT compositions in long format for ggplot
rawComp4 <- pivot_longer(rawComp,names_to="nt",values_to="composition",c(A,T,G,C,N)) %>%
	mutate(nt=factor(nt,levels=c("A","G","C","T","N"))) #set order for nucleotide cols
rawComp2 <- pivot_longer(rawComp,names_to="nt",values_to="composition",c(AT,GC,N))

# Process reverse complement compositions for right zoomed-in plots
if(verbose) cat(TimeDiff(time0),"Executing pre-trim seqtk reverse pass ... ")
timeStep <- Sys.time()
seqtkBuffer <- ExecSeqtk(seqFile,revComp=TRUE,head=headFastq) #compute compositions via seqtk
if(verbose) cat("done",TimeDiff(timeStep),"\n")
revComp <- ExtractCompositions(seqtkBuffer,percHide) #extract revComp fastq compositions from seqtkBuffer
rm(seqtkBuffer) #destroy seqtk buffer for memory
# Transform revComp NT compositions in long format for ggplot
revComp4 <- pivot_longer(revComp,names_to="nt",values_to="composition",c(A,T,G,C,N)) %>%
	mutate(nt=factor(nt,levels=c("A","G","C","T","N"))) #set order for nucleotide cols
revComp2 <- pivot_longer(revComp,names_to="nt",values_to="composition",c(AT,GC,N))


## Generate graph objects for 1st round ##

# Process full compositions for center plots. Pushed to --noTrim condition.
# areaComp4 <- AreaPlot(rawComp4)
# areaComp2 <- AreaPlot(rawComp2)
# Process 5' end compositions for left zoomed-in plots
areaComp4Zoom5 <- AreaPlot(rawComp4,zoom5Len)
areaComp2Zoom5 <- AreaPlot(rawComp2,zoom5Len)
# Process 3' end with reverse compositions for right zoomed-in plots
areaComp4Zoom3 <- AreaPlot(revComp4,zoom3Len)
areaComp2Zoom3 <- AreaPlot(revComp2,zoom3Len)


## Open connection to a markdown file
mdFile <- paste0(tmpDir,"/",outPrefix,".md")
mdReport <- paste0(workDir,"/",outPrefix,".html")
mdConn <- file(mdFile,open="wt")
writeLines(paste0('---\ntitle: "SNIKT contamination report"\n---\n'),mdConn) #make metadata section, add title
writeLines(paste0("## Input: ",seqFileName,"\n\n---\n\n"),mdConn)

# RenderMd()


## Condition path when user asks to not trim ##
# Program quits after following this block
if(!is.na(trim5Len) & !is.na(trim3Len) & trim5Len==0 & trim3Len==0){
	if(verbose) cat(TimeDiff(time0),"Preparing graphs without trimming ... ")
	timeStep <- Sys.time()
	# Process full compositions for center plots
	areaComp4 <- AreaPlot(rawComp4,xBreaks=xBreaks)
	areaComp2 <- AreaPlot(rawComp2,xBreaks=xBreaks)
	areaLen <- dim(rawComp4)[1]/5
	noTrimGrid <- CreateGrid(legends,areaComp4Zoom5,areaComp4Zoom3,areaComp2Zoom5,areaComp2Zoom3,areaComp4,areaComp2,areaLen)
	noTrimFile <- paste0(tmpDir,"/",outPrefix,"-noTrim.png")
	ggsave(noTrimFile,plot=noTrimGrid,height=200,width=400,units="mm")
	if(verbose) cat("done",TimeDiff(timeStep),"\n")
	# Build a section title depending on whether the fastq was subsetted
	sectionTitle <- paste0("### Nucleotide compositions without trimming")
	if(headFastq>0){sectionTitle <- paste0(sectionTitle," - up to top ",headFastq," reads")}
	sectionTitle <- paste0(sectionTitle,"\n\n")
	writeLines(sectionTitle,mdConn)
	# writeLines(paste0("Fastq: ",seqFileName,"\n\n"),mdConn)
	writeLines(knitr::kable(summRawFq),mdConn) #write raw fastq summary to md
	writeLines(paste0("\n\n![](",noTrimFile,")\n\n"),mdConn)
	# generate report with a function
	# quit
	RenderMd(mdConn,mdFile,mdReport)
}


## Default condition path with trimming ##


## Interactive mode - take user input
if(is.na(trim5Len)|is.na(trim3Len)){
	# Generate temporary graph and save
	userGrid4 <- CreateGrid(legends,areaComp4Zoom5,areaComp4Zoom3,areaComp2Zoom5,areaComp2Zoom3)
	interGraphFile <- paste0(workDir,"/",outPrefix,"-",TimeStamp(),".png")
	ggsave(interGraphFile,plot=userGrid4,height=200,width=400,units="mm")
	# Send messages to user
	cat("SNIKT - interactive mode")
	cat("\nThe temporary graphs are ready in:",interGraphFile)
	cat("\nPlease check the graphs to decide on trim lengths.")
	cat("\n\nIf you are unsatisfied with the zoom window, please re-run with custom -Z and -z values. To quit, hit Enter after feeding SIGTERM (Ctrl+c).")
	cat(paste0("\n\nFilter length is currently set to drop any post-trim reads below length: ",filterLen))
	cat("\nIf you want to decide appropriate filter length, please quit (Ctrl+c,Enter) and rerun with -n to get a 6-panel graph that helps decide both trim and filter values.")
	cat("\n\nWhen you are ready, please fill-in below.\n")
	# Read user input
	# trim5Max <- as.integer(summRawFq$`Max Length`)
	trim5Max <- zoom5Len #set max possible 5' trim to zoom5 window size for interactive mode
	trim5Len <- ReadTrimLen(paste0("Your choice for 5' trimming length [0-",trim5Max,"]: "),trim5Max)
	# Trim 3 max is the smaller of what's left after trim 5 from max or zoom 3 size
	trim3Max <- as.integer(summRawFq$`Max Length`)-trim5Len
	if(trim3Max>zoom3Len) trim3Max <- zoom3Len
	# trim3Max <- zoom3Len #setting the max possible value to the zoom 3 window in interactive mode
	trim3Len<-ReadTrimLen(paste0("Your choice for 3' trimming length [0-",trim3Max,"]: "),trim3Max)
	if(trim5Len+trim3Len>as.integer(summRawFq$`Max Length`)){stop("ERR: Sum of trimming lengths is larger than max read length.")} #verification that sum of trim lengths < Max read length
	
	# Delete intermediate graph file if not needed anymore
	if(!file.remove(interGraphFile)){write(paste0("WARN: Could not remove temporary graph file: ",interGraphFile),stderr())}
}


## Make 1st set of graphs with trim lines
areaComp4Zoom5 <- areaComp4Zoom5+geom_vline(xintercept=trim5Len,linetype="dashed")
areaComp2Zoom5 <- areaComp2Zoom5+geom_vline(xintercept=trim5Len,linetype="dashed")
areaComp4Zoom3 <- areaComp4Zoom3+geom_vline(xintercept=trim3Len,linetype="dashed")
areaComp2Zoom3 <- areaComp2Zoom3+geom_vline(xintercept=trim3Len,linetype="dashed")
if(verbose) cat(TimeDiff(time0),"Preparing graphs pre-trimming ... ")
timeStep <- Sys.time()
preTrimGrid <- CreateGrid(legends,areaComp4Zoom5,areaComp4Zoom3,areaComp2Zoom5,areaComp2Zoom3)
preTrimFile <- paste0(outPrefix,"-preTrim.png")
ggsave(preTrimFile,plot=preTrimGrid,height=200,width=400,units="mm")
if(verbose) cat("done",TimeDiff(timeStep),"\n")

# Export pre-trim graphs in markdown
# Build a section title depending on whether the fastq was subsetted
sectionTitle <- paste0("### Nucleotide compositions pre-trimming")
if(headFastq>0){sectionTitle <- paste0(sectionTitle," - up to top ",headFastq," reads")}
sectionTitle <- paste0(sectionTitle,"\n\n")
writeLines(sectionTitle,mdConn)
# writeLines(paste0("Fastq: ",seqFileName,"\n\n"),mdConn)
writeLines(knitr::kable(summRawFq),mdConn) #write raw fastq summary to md
writeLines(paste0("\n\n![](",preTrimFile,")\n\n"),mdConn)
if(filterLen>0){writeLines(paste0("Filter: dropped reads shorter than ",filterLen," bases.\n\n"),mdConn)}


## Trim reads with seqtk
# if(threads<=1){
# 	trimFile <- TrimFiltSeq(seqFile,trim5Len,trim3Len,filterLen) #serial trimming
# } else{
# 	# This will need decompression handling.
# 	trimFile <- TrimSeqParallel(seqFile,trim5Len,trim3Len,threads) #parallel trimming
# }
trimFile <- paste0(workDir,"/",trimFileName)
trimFile <- TrimFiltSeq(seqFile,trim5Len,trim3Len,filterLen,trimFile)
cat("SNIKT cleaned reads are available in: ",trimFile,"\n")


## Execute 2nd round of seqtk and compute compositions for trimmed reads
# Process full-fastq base compositions for full-center and zoomed-left plots
if(verbose) cat(TimeDiff(time0),"Executing post-trim seqtk forward pass ... ")
timeStep <- Sys.time()
seqtkBuffer <- ExecSeqtk(trimFile) #compute compositions via seqtk
if(verbose) cat("done",TimeDiff(timeStep),"\n")

# trimFileName <- basename(trimFile)
# summTrimFq <- ExtractSummary(seqtkBuffer,trimFileName) #extract fastq summary from seqtkBuffer
summTrimFq <- ExtractSummary(seqtkBuffer) #extract fastq summary from seqtkBuffer
rawComp <- ExtractCompositions(seqtkBuffer,percHide)  #extract fastq compositions from seqtkBuffer
# Transform NT compositions in long format for ggplot
rawComp4 <- pivot_longer(rawComp,names_to="nt",values_to="composition",c(A,T,G,C,N)) %>%
	mutate(nt=factor(nt,levels=c("A","G","C","T","N"))) #set order for nucleotide cols
rawComp2 <- pivot_longer(rawComp,names_to="nt",values_to="composition",c(AT,GC,N))

# Process full-fastq reverse complement compositions for zommed-right plots
if(verbose) cat(TimeDiff(time0),"Executing post-trim seqtk reverse pass ... ")
timeStep <- Sys.time()
seqtkBuffer <- ExecSeqtk(trimFile,revComp=TRUE) #compute compositions via seqtk
if(verbose) cat("done",TimeDiff(timeStep),"\n")
revComp <- ExtractCompositions(seqtkBuffer,percHide)  #extract revComp fastq compositions from seqtkBuffer
rm(seqtkBuffer) #destroy seqtk buffer for memory
# Transform revComp NT compositions in long format for ggplot
revComp4 <- pivot_longer(revComp,names_to="nt",values_to="composition",c(A,T,G,C,N)) %>%
	mutate(nt=factor(nt,levels=c("A","G","C","T","N"))) #set order for nucleotide cols
revComp2 <- pivot_longer(revComp,names_to="nt",values_to="composition",c(AT,GC,N))


## Make 2nd set of graphs from trimmed read compositions
if(verbose) cat(TimeDiff(time0),"Preparing graphs post-trimming ... ")
timeStep <- Sys.time()
# Process full compositions for center plots.
areaComp4 <- AreaPlot(rawComp4,xBreaks=xBreaks)
areaComp2 <- AreaPlot(rawComp2,xBreaks=xBreaks)
# Process 5' end compositions for left zoomed-in plots
areaComp4Zoom5 <- AreaPlot(rawComp4,zoom5Len)
areaComp2Zoom5 <- AreaPlot(rawComp2,zoom5Len)
# Process 3' end with reverse compositions for right zoomed-in plots
areaComp4Zoom3 <- AreaPlot(revComp4,zoom3Len)
areaComp2Zoom3 <- AreaPlot(revComp2,zoom3Len)
areaLen <- dim(rawComp4)[1]/5
postTrimGrid <- CreateGrid(legends,areaComp4Zoom5,areaComp4Zoom3,areaComp2Zoom5,areaComp2Zoom3,areaComp4,areaComp2,areaLen)
postTrimGraphFile <- paste0(tmpDir,"/",outPrefix,"-noTrim.png")
ggsave(postTrimGraphFile,plot=postTrimGrid,height=200,width=400,units="mm")
if(verbose) cat("done",TimeDiff(timeStep),"\n")

# Export post-trim graphs in markdown
writeLines(paste0("<br>\n\n---\n\n<br>\n\n### Nucleotide compositions post-trimming\n\n"),mdConn)
# writeLines(paste0("---\n\n---\n\n### Nucleotide compositions post-trimming\n\n"),mdConn)
writeLines(paste0("**Fastq: ",basename(trimFile),"**\n\n"),mdConn)
writeLines(knitr::kable(summTrimFq),mdConn) #write raw fastq summary to md
writeLines(paste0("\n\n![](",postTrimGraphFile,")\n\n"),mdConn)

# Render report and finish
RenderMd()
