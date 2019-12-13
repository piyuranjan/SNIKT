#!/usr/bin/env Rscript

## Commenting suggestions ##
## Section headers like this ##
## Content headers like this
### Full-line comment like this
### and can continue on second line like this
###   or this
#x<-as.numeric(y); #inline comment like this

startTime<-Sys.time()
packages<-c("tidyverse","grid","gridExtra","argparser","lubridate")
tmp<-lapply(packages,function(x) suppressPackageStartupMessages(require(x,character.only=T)))
# options(readr.num_columns=0) #suppress runtime output for read_delim and other functions

ArgParser2<-function(description,name=NULL) 
	{
	### This function is borrowed from the arg_parser() function of the argparser library
	### that creates an arg_parser object with modifications. For original functionality check ?arg_parser
	### options -- and --opts are commented in this section for their exclusion from the help document
	if(is.null(name))
		{
		prefix<-"--file="
		name<-sub(prefix,"",grep(paste(prefix,"(.+)",sep=""),commandArgs(),value=TRUE))
		}
	if(length(name)==0){name<-"<script>"}
	parser<-structure(list(name=name,description=description),class="arg.parser")
	# parser<-add_argument(parser,"--","placeholder",flag=TRUE)
	parser<-add_argument(parser,"--help","show this help message and exit",flag=TRUE)
	# parser<-add_argument(parser,"--opts","RDS file containing argument values",short="-x")
	parser
	}

Preprocess<-function(compFile,percentFilter=0.005) #This function is deprecated and will be removed in a future commit
	{
	### This function stores nucleotide composition per base data from a file in a data frame
	### subsets the data frame by a user given threshold to remove long tail, which is often extremely noisy
	### adds AT, GC columns to it
	# maxBases<-read_delim(compFile,"\t",skip=3,n_max=1)[,2]
	maxBases<-as.numeric(read_delim(compFile,"\t",skip=1,n_max=1)[,2])
	numFilter<-maxBases*percentFilter #number of base positions to subset until, avoids rendering noise at the end
	# df<-read_delim(compFile,"\t",skip=3)
	# colnames(df)<-c("POS","bases","A","C","G","T","N","avgQ","errQ","low","high")
	df<-read_delim(compFile,"\t")
	df<-dplyr::filter(df,bases>=numFilter) #remove low confidence compositions by filtering the long tail
	df<-mutate(df,AT=A+T,GC=G+C)
	}

ComputeComp<-function(seqFile,percentFilter=0.005,revComp=F)
	{
	### This function computes nucleotide composition per base using the program seqtk
	### stores the compositions in a data frame
	### subsets the data frame by a user given threshold to remove long tail, which is often extremely noisy
	### adds AT, GC columns to it
	
	## Prepare command and execute seqtk
	if(revComp==F) #compute forward 5' aligned compositions
		{command<-paste("seqtk fqchk",seqFile)}
	else #compute reverse complement for 3' alignment and calculate compositions
		{command<-paste("seqtk seq -r",seqFile,"|seqtk fqchk -")}
	compositions<-system(command,intern=T) #execute seqtk on shell and retrieve data
	
	## Process the seqtk output, store and subset
	colNames<-c("POS","bases","A","C","G","T","N","avgQ","errQ","low","high") #full list of columns in seqtk output
	df<-read_delim(compositions,"\t",skip=3,col_names=colNames)
	df<-df[,1:8] #leaving out last 3 unneccesary columns
	
	## Filter compositions by percentage and add paired compositions
	maxBases<-as.numeric(df[1,2])
	numFilter<-maxBases*percentFilter #number of base positions to subset until, avoids rendering noise at the end
	df<-dplyr::filter(df,bases>=numFilter) #remove low confidence compositions by filtering the long tail
	df<-mutate(df,AT=A+T,GC=G+C) #add paired AT and GC compositions
	}

ComputeTrim<-function(df,zoomLen)
	{
	### This function is a work in progress. Determine an algorithm to compute the Trim threshold.
	trim<-max(which(df[0:(zoomLen/2), "avgQ"] < mean(df[300:500,]$avgQ))) + 1
	return(trim)
	}

format_si <- function(...) {
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
    prefix <- c("y", "z", "a",  "f", "p", "n", "µ",  "m",   "",  "K", "M",  "G",  "T",  "P",  "E", "Z", "Y")
    
    # Vector with array indices according to position in intervals
    i<-findInterval(abs(x),limits)
    
    # Set prefix to " " for very small values < 1e-24
    i<-ifelse(i==0,which(limits==1e0),i)
    paste(format(round(x/limits[i],1),trim=TRUE,scientific=FALSE,...),prefix[i])
  }
}

AreaPlot<-function(df,zoomLen=NULL,trim=NULL,pad=0.025)
	{
	### This function makes an area plot for given compositions per position
	### adds average base quality and sequence length distribution for full length plots
	### subsets data frame for zoomed plots
	
	if(!is.null(zoomLen)){df<-dplyr::filter(df,POS<=zoomLen)} #subsetting condition to follow for zoom plots
	
	ap<-ggplot(df,aes(x=POS,y=composition,fill=nt))+
		geom_area(alpha=0.75)+ #area plot; lowered alpha to see gridlines
		geom_line(aes(y=avgQ,color="Average\nPhred\nScore"),alpha=0.75)+ #avg Phred score
	  geom_line(aes(y=bases/max(bases)*100,color="Nucleotide\nPer\nPosition"))+ #make normalized length freq as line
	  scale_y_continuous(labels=paste0(c(0,25,50,75,100),"%"),expand=c(pad,0), #add % label, adjust top/bot plot padding
      sec.axis=sec_axis(~.*max(df$bases)/100,name="Frequency",labels=format_si()))+ #secondary axis uses SI scale
		theme_bw()+ #remove grey background
	  theme(title=element_blank(),axis.line=element_line(color="black"),axis.line.y.right=element_line(color="red"), 
      axis.ticks.y.right=element_line(color="red"),axis.text.y.right=element_text(color="red")) #change sec axis line, ticks, text to red
	
	if(is.null(zoomLen)) #condition to skip for zoom plots ## remove commented code in future version
		{
		ap<-ap+#geom_line(aes(y=bases/max(bases)*100,color="Nucleotide\nPer\nPosition"))+ #make normalized length freq as line
			# scale_y_continuous(labels=paste0(c(0,25,50,75,100),"%"),expand=c(pad,0), #add % label, adjust top/bot plot padding
			                   # sec.axis=sec_axis(~.*max(df$bases)/100,name="Frequency"))+ #rev-norm freq for axis
		  scale_x_continuous(expand=c(pad,0))+ # adjust left/right plot padding
		  # theme(axis.line=element_line(color="black"),axis.line.y.right=element_line(color="red"), 
		    # axis.ticks.y.right=element_line(color="red"),axis.text.y.right=element_text(color="red"))+ #changes secondary (right) axis line, ticks, text to red
		  scale_color_manual(name=element_blank(),breaks=c("Average\nPhred\nScore","Nucleotide\nPer\nPosition"),values=c("Average\nPhred\nScore"="black","Nucleotide\nPer\nPosition"="red"))+ #create line legend
		  guides(fill=guide_legend(title=NULL,order=1),color=guide_legend(order=2)) #enforce legends appearing in consistent order
		return(ap) #return plot with legends for extraction; legends and axis labels will be removed prior to assembling plot grid
		}
	
	### positioned here to avoid null zoomLen comparison; from zoom argument, determine if forward or reverse scale applies
	if(deparse(substitute(zoomLen))=="zoom5Len") {
	  ap<-ap+scale_x_continuous(expand=c(pad,0)) # adjust left/right plot padding
	} else {
	  ap<-ap+scale_x_reverse(expand=c(pad,0)) # reverse x-axis, adjust left/right plot padding
	}
	
	ap<-ap+scale_color_manual(values=c("Average\nPhred\nScore"="black","Nucleotide\nPer\nPosition"="red"))+ #fix avqQ to black
	  # scale_y_continuous(labels=paste0(c(0,25,50,75,100),"%"),expand=c(pad,0))+ # add % label, adjust top/bottom plot padding
	  theme(legend.position="none") #remove all title elements, legend from zoom plot. add left/bot axis lines
	
	if(!is.null(trim)){ap<-ap+geom_vline(xintercept=trim,linetype="dashed")} #add dashed vertical line at trim position
	
	return(ap) #return zoom plot
	}

##########
## Main ##
##########

### Status message that records time for libraries to load. Feel free to uncomment the following line if you want to see it.
# cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Libraries loaded\n",sep="")

## User defined parameter parsing ##

## Defining and reading in user arguments
args<-ArgParser2("FastQ QC and sequence over-representation check. \nCalculate and plot per-position nucleotide composition, for finding \nadapter contamination/over-representation in sequences. \nAuthors: Piyush Ranjan, Christopher Brown\n\nFor more info, please check: https://github.com/piyuranjan/MetagenomeAdapterIdentify") #defines the arg_parser object
args<-add_argument(args,"fastq","FastQ file for composition analysis; required",default=NULL)
args<-add_argument(args,"--noTrim","Disable positional trimming; useful for short-read data",flag=T,type="logical",short="-n")
args<-add_argument(args,"--trim5Len","[1..(maxSeqLen-trim3Len)] NT to trim from aligned 5' side; interactive mode if unset",type="numeric",short="-t5")
args<-add_argument(args,"--trim3Len","[1..(maxSeqLen-trim5Len)] NT to trim from aligned 3' side; interactive mode if unset",type="numeric",short="-t3")
args<-add_argument(args,"--filter","[0..1] Filter the composition tail by a fraction",default=0.005,short="-f")
args<-add_argument(args,"--bzoom","[1..maxSeqLen] NT to zoom in from aligned 5' begin",default=300,short="-b")
args<-add_argument(args,"--ezoom","[1..maxSeqLen] NT to zoom in from aligned 3' end",default=100,short="-e")
args<-add_argument(args,"--outfile","[file.ext] File for saving graphs [default: fastqNoExt-fastqCompositions.png]",short="-o")
args<-add_argument(args,"--verbose","Enable status messages",flag=T,type="logical",short="-v")
args<-parse_args(args,argv=commandArgs(trailingOnly=T))
if(!file.exists(args$fastq)){stop("Need a sequence file to proceed. See help with -h.")} #condition to kill if no seqFile

### Resolving bug regarding generation of Rplots.pdf
if(file.exists("Rplots.pdf")) #renaming any preexisting Rplots.pdf so that it doesn't overwrite
	{invisible(file.rename("Rplots.pdf",paste0("Rplots_",str_replace_all(Sys.time(),"[-:]|\\h+",""),".pdf")))}

## Associating user flags with variables
seqFile<-args$fastq
noTrim<-args$noTrim
trim5Len<-args$trim5Len
trim3Len<-args$trim3Len
percentFilter<-args$filter #percent bases to filter until; parameter that avoids rendering noise at the end
zoom5Len<-args$bzoom
zoom3Len<-args$ezoom
verbose<-args$verbose
outFile<-sub("\\.f[ast]{0,3}(a|q)(\\.gz)?","-fast\\1Compositions.png",seqFile,perl=T) #file to export graphs in
if(!is.na(args$outfile)){outFile<-args$outfile} #override filename if supplied by user
if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Libraries loaded and user parameters set\n",sep="")}

## Compute nucleotide compositions, quality estimates and generate plots ##

## Process regular full base compositions for center plots
# rawComp<-Preprocess(rawQualFile,percentFilter) #sending to function for file reading
rawComp<-ComputeComp(seqFile,percentFilter) #compute compositions via seqtk
rawComp4nt<-gather(rawComp,key=nt,value=composition,c(A,T,G,C,N)) #gathers, reorders nt in long format to use with ggplot
rawComp4nt$nt<-factor(rawComp4nt$nt,levels=c("A","G","C","T","N")) #setting the order for nucleotide cols
areaRawComp4nt<-AreaPlot(rawComp4nt) #plotting full area graph with 4 nt
rawComp2nt<-gather(rawComp,key=nt,value=composition,c(AT,GC,N)) #reorders AT, GC in long format for ggplot
areaRawComp2nt<-AreaPlot(rawComp2nt) #plotting full area graph with combination compositions

## Extract legends, then remove from plots
nt4Grob<-ggplotGrob(areaRawComp4nt)$grobs
nt2Grob<-ggplotGrob(areaRawComp2nt)$grobs
nt4Leg<-nt4Grob[[which(map_chr(nt4Grob,function(x) x$name)=="guide-box")]]
nt2Leg<-nt2Grob[[which(map_chr(nt2Grob,function(x) x$name)=="guide-box")]]
areaRawComp4nt<-areaRawComp4nt+theme(legend.position="none")
areaRawComp2nt<-areaRawComp2nt+theme(legend.position="none")
if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Forward full compositions calculated and graphed\n",sep="")}

## Process 5' end compositions for left zoomed-in plots
areaRawComp4ntZoom5<-AreaPlot(rawComp4nt,zoom5Len,ComputeTrim(rawComp4nt,zoom5Len))
areaRawComp2ntZoom5<-AreaPlot(rawComp2nt,zoom5Len,ComputeTrim(rawComp2nt,zoom5Len))
if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Forward zoomed compositions calculated and graphed\n",sep="")}

## Process reverse complement compositions for right zoomed-in plots
# revComp<-Preprocess(revQualFile,percentFilter) #reading and filtering reverse complement compositions
revComp<-ComputeComp(seqFile,percentFilter,revComp=T)
revComp4nt<-gather(revComp,key=nt,value=composition,c(A,T,G,C,N)) #gathers, reorders nt in long format to use with ggplot
revComp4nt$nt<-factor(revComp4nt$nt,levels=c("A","G","C","T","N")) #setting the order for nucleotide cols
areaRevComp4ntZoom3<-AreaPlot(revComp4nt,zoom3Len,ComputeTrim(revComp4nt,zoom3Len))
revComp2nt<-gather(revComp,key=nt,value=composition,c(AT,GC,N)) #reorders AT, GC in long format for ggplot
areaRevComp2ntZoom3<-AreaPlot(revComp2nt,zoom3Len,ComputeTrim(revComp2nt,zoom3Len))
if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Reverse zoomed compositons calculated and graphed\n",sep="")}

## Export graphs in figure ##

## Setting graphs in grids
leftPlots<-arrangeGrob(areaRawComp4ntZoom5,areaRawComp2ntZoom5,ncol=1,top=textGrob("Aligned 5' beginning",gp=gpar(fontface=2,fontsize=14)))
centerPlots<-arrangeGrob(areaRawComp4nt,areaRawComp2nt,ncol=1,
						top=textGrob(paste("Compositions to length",dim(rawComp)[1],"aligned 5'"),gp=gpar(fontface=2,fontsize=14)),
						right=textGrob(label="Read frequency (NT per position)",rot=90,gp=gpar(fontface=2,fontsize=14,col="red")))
rightPlots<-arrangeGrob(areaRevComp4ntZoom3,areaRevComp2ntZoom3,ncol=1,
						top=textGrob("Aligned 3' ending",gp=gpar(fontface=2,fontsize=14)))
areaArranged<-arrangeGrob(leftPlots,centerPlots,rightPlots,ncol=3,widths=c(1,2,1),
							top=textGrob(label="Sequence Representation Per Position",gp=gpar(fontface=2,fontsize=20)),
							left=textGrob(label="Nucleotide composition (0-100), Average Phred score",rot=90,gp=gpar(fontface=2,fontsize=14)),
							bottom=textGrob(label="Read length (NT positions)",gp=gpar(fontface=2,fontsize=14)))
finalPlotGrid<-grid.arrange(areaArranged,arrangeGrob(nt4Leg,nt2Leg),ncol=2,widths=c(28,2))
if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] All graphs prepared in grid\n",sep="")}

## Export graphs to file
ggsave(outFile,plot=finalPlotGrid,height=200,width=400,units="mm")
### Remove extra Rplots.pdf if generated by this instance; relates to a bug
if(file.exists("Rplots.pdf")){invisible(file.remove("Rplots.pdf"))}

if(isTRUE(verbose)){cat("[",round(time_length(Sys.time()-startTime,unit="second"),0),"s] Graph grid written to file\n",sep="")}
