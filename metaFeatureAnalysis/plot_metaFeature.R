# Description: This function is for plotting metaorf
#
# Author:      Mohamood Adhil
# Date:        24th May 2017 
# help:        Rscript <R-script> -h (or) --help
# TODO:        
# #################################################################################
# 
.libPaths("/storage/home/amohamme/Rlibrary")
library(argparse)
##############

parser <- ArgumentParser()
parser <- ArgumentParser(description='This script is for plotting metaORF values')
parser$add_argument("-mf", dest="inmet", help="Comma seperated text files (Mandatory)", required = TRUE)
parser$add_argument("-tg", dest="tags", help="semi colan (;) seperated within the metaorf txt file and comma (,) seperated between the metaorf files (Example:SampleASubtype1;SampleASubtype2,SampleBSubtype1;SampleBSubtype2), If not provided all the subtypes will be used",  default="None")
parser$add_argument("-rtype", dest="rtype", help="ratio type l2fc=log2(IP/Input) or just ratio=(IP/Input), if l2fc then exponential function is used to convert all values to positive scale (Default = l2fc)", default="l2fc", choices=c("l2fc","ratio"))
parser$add_argument("-sm", dest="sm", help="smoothing the curve using geom_smooth function (Default=yes)", default="yes",choices=c("yes","no"))
parser$add_argument("-cf", dest="cf", help="95 percent confidence interval to be plotted (Default=no)", default="no",choices=c("yes","no"))
parser$add_argument("-nb", dest="nbases", help="Number of bases used for merging (Mandatory)", required = TRUE)
parser$add_argument("-sf", dest="stag", help="suffix tag for plot (Mandatory)", required = TRUE)
parser$add_argument("-ma", dest="main", help="Main Title (Mandatory)", required = TRUE)
parser$add_argument("-dir", dest="outdir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

####################################################################################

#Library
library(ggplot2)
library(gridExtra)
library(plyr)
set.seed(1)

#Arguments
inmet = args$inmet
metfiles = strsplit(inmet,",")[[1]]

tags = args$tags
if (tags != "None"){
smtags = strsplit(tags,",")[[1]]
}

rtype = args$rtype

smooth = args$sm

cfi = args$cf

nbases = as.numeric(args$nbases)

stag = args$stag

main = args$main

outdir = args$outdir
dir.create(outdir, showWarnings = FALSE)

##
metalist = list()
cli = list()
pli = list()
rname = c()
cn = 0
allt = c()
for (i in 1:length(metfiles)){
  alldata = read.csv(metfiles[i],sep="\t",header=TRUE)
  print (metfiles[i])
  if (tags == "None"){
  tago = unique(as.vector(alldata$subtype))
  } else {
  tago = strsplit(smtags[i],";")[[1]]
  }
  for (j in 1:length(tago)){
  cn = cn+1
  dat = alldata[which(alldata$subtype == tago[j]),]
  allt = c(allt,tago[j])
  metalist[[cn]] = dat
  }
}
metaorfall = ldply(metalist, data.frame)
metaorf = metaorfall[,c(1,2,4,5,6)]
colnames(metaorf) = c("Distance","likelihood","Sample","CF_low","CF_high")
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
metaorf$likelihood = as.numeric(metaorf$likelihood)
metaorf$CF_low = as.numeric(metaorf$CF_low)
metaorf$CF_high = as.numeric(metaorf$CF_high)

if (rtype == "l2fc"){
metaorf$likelihood = exp(as.numeric(metaorf$likelihood))
metaorf$CF_low = exp(as.numeric(metaorf$CF_low))
metaorf$CF_high = exp(as.numeric(metaorf$CF_high))
}

if (smooth == "yes"){
if (cfi == "yes"){
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + ggtitle(main) +
 geom_ribbon(aes(ymin=metaorf$CF_low, ymax=metaorf$CF_high), alpha=0.1) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Average Score") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

 fname = paste(outdir,"/MetaFeature_Score_",stag,"_no.jpeg",sep="")
 ggsave(p, file=fname, width=10, height=6, dpi=1000)
}
if (cfi == "no"){
 p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Average Score") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

 fname = paste(outdir,"/MetaFeature_Score_",stag,"_no.jpeg",sep="")
 ggsave(p, file=fname, width=10, height=6, dpi=1000)
}
metaorf = metaorfall[,c(1,3,4)]
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Average Count") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

fname = paste(outdir,"/MetaFeature_Count_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)
}

if (smooth == "no"){
if (cfi == "yes"){
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_line() +
  geom_ribbon(aes(ymin=metaorf$CF_low, ymax=metaorf$CF_high), alpha=0.1) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Normalized bTMP Score (IP/Input)") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))
  fname = paste(outdir,"/MetaFeature_Score_",stag,"_no.jpeg",sep="")
  ggsave(p, file=fname, width=10, height=6, dpi=1000)
}
if (cfi == "no"){
 p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_line() +
  ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Normalized bTMP Score (IP/Input)") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))
  fname = paste(outdir,"/MetaFeature_Score_",stag,"_no.jpeg",sep="")
  ggsave(p, file=fname, width=10, height=6, dpi=1000)
}
metaorf = metaorfall[,c(1,3,4)]
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_line()  + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Average Count") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

fname = paste(outdir,"/MetaFeature_Count_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)
}
