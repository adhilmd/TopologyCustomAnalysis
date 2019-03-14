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
parser$add_argument("-mf", dest="inmet", help="Comma seperated metaorf data files (RData file) (Mandatory)", required = TRUE)
parser$add_argument("-tg", dest="tags", help="semi colan (;) seperated within the metaorf RData file and comma (,) seperated between the metaorf files (Mandatory)", required = TRUE)
parser$add_argument("-nt", dest="newt", help="semi colan (;) seperated within the metaorf RData file and comma (,) seperated between the metaorf files (Mandatory)", required = TRUE)
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
smtags = strsplit(tags,",")[[1]]

newt = args$newt
nwtags = strsplit(newt,",")[[1]]

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
  print (metfiles[i])
  load(metfiles[i])
  alldata = finallist[[2]]
  tago = strsplit(smtags[i],";")[[1]]
  tagr = strsplit(nwtags[i],";")[[1]]
  for (j in 1:length(tago)){
  print (tago[j])
  cn = cn+1
  dat = alldata[which(alldata$sample == tago[j]),]
  dat$sample = tagr[j]
  allt = c(allt,tagr[j])
  metalist[[cn]] = dat
  }
}
metaorf = ldply(metalist, data.frame)
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$likelihood = exp(as.numeric(metaorf$likelihood))
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + theme_grey(base_size = 20) + ggtitle(main) + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14, family = "arial")) + 
  theme(legend.title = element_text(size=10, face="bold", family = "arial"))  + theme(legend.text = element_text(size = 8, family = "arial")) + 
  xlab("Average Distance") + ylab("Average Ratio IP/Input") + geom_vline(xintercept = 0) + geom_vline(xintercept = 1000) + 
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12, family="arial"))

fname = paste(outdir,"/MetaORF_Intensity_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)

metalist = list()
cli = list()
pli = list()
rname = c()
cn = 0
allt = c()
for (i in 1:length(metfiles)){
  print (metfiles[i])
  load(metfiles[i])
  alldata = finallist[[3]]
  tago = strsplit(smtags[i],";")[[1]]
  tagr = strsplit(nwtags[i],";")[[1]]
  for (j in 1:length(tago)){
  print (tago[j])
  cn = cn+1
  dat = alldata[which(alldata$sample == tago[j]),]
  dat$sample = tagr[j]
  allt = c(allt,tagr[j])
  metalist[[cn]] = dat
  }
}
metaorf = ldply(metalist, data.frame)
colnames(metaorf) = c("Distance","likelihood","Sample")
metaorf$Sample = factor(metaorf$Sample, levels=c(allt))
p <- ggplot(data= metaorf, aes(x=Distance, y=likelihood, colour=Sample)) + geom_smooth(se=F) + theme_grey(base_size = 20) + ggtitle(main) +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=14, family = "arial")) +
  theme(legend.title = element_text(size=10, face="bold", family = "arial"))  + theme(legend.text = element_text(size = 8, family = "arial")) +
  xlab("Average Distance") + ylab("Average Gene Density") + geom_vline(xintercept = 0) + geom_vline(xintercept = 1000) +
  geom_hline(yintercept=0) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title = element_text(face="bold", size = 12, family="arial"))

fname = paste(outdir,"/MetaORF_Count_",stag,"_no.jpeg",sep="")
ggsave(p, file=fname, width=10, height=6, dpi=1000)
