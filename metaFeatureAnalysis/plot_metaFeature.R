# Description: This function is for plotting metaorf
#
# Author:      Mohamood Adhil
# Date:        24th May 2017 
# help:        Rscript <R-script> -h (or) --help
# TODO:        
# #################################################################################
# 
library(argparse)
##############

parser <- ArgumentParser()
parser <- ArgumentParser(description='This script is for plotting metaORF values')
parser$add_argument("-mf", dest="inmet", help="Comma seperated text files (Mandatory)", required = TRUE)
parser$add_argument("-tg", dest="tags", help="semi colan (;) seperated within the metaorf txt file and comma (,) seperated between the metaorf files (Example:SampleASubtype1;SampleASubtype2,SampleBSubtype1;SampleBSubtype2), If not provided all the subtypes will be used",  default="None")
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
metaorf = metaorfall[,c(1,2,5)]
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$likelihood = exp(as.numeric(metaorf$likelihood))
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) + 
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) + 
  xlab("Average Distance") + ylab("Average Score") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) + 
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

fname = paste(outdir,"/MetaFeature_Score_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)

metaorf = metaorfall[,c(1,3,5)]
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14)) +
  theme(legend.title = element_text(size=10, face="bold"))  + theme(legend.text = element_text(size = 8)) +
  xlab("Average Distance") + ylab("Average Feature Count") + geom_vline(xintercept = 0) + geom_vline(xintercept = nbases) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12))

fname = paste(outdir,"/MetaFeature_Count_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)
