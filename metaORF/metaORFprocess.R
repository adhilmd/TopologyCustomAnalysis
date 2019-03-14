# Description: This function is to average metaORF values for a given length and create RData file
#
# Author:      Mohamood Adhil
# Date:        20th April 2017
# Usage:       Rscript metaORFprocess.R 
# Example1:    Rrun 
# help:        Rscript <R-script> -h (or) --help
# TODO:        
# #################################################################################

#library
library(argparse)
##############

parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to average metaORF values for a given length and create RData file')
parser$add_argument("-i", dest="inpsamples", help="Input directory where all the bedgraph files are present (Mandatory)", required = TRUE)
parser$add_argument("-sm", dest="samplenames", help="comma seperated sample names", required = TRUE)
parser$add_argument("-st", dest="st", help="strand yes or no", required = TRUE)
parser$add_argument("-ft", dest="ft", type="integer", help="Average number of bases for feature Default = 1000", default=1000)
parser$add_argument("-fb", dest="flankingbases", type="integer", help="Flanking number of bases for feature Default = 500", default=500)
parser$add_argument("-pr", dest="prefix", required=TRUE, help="Prefix tag for output (RData) file (Mandatory)")
parser$add_argument("-dir", dest="outdir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

###################################################################
library(plyr)
library(dplyr)
library(tidyr)
library(zoo)

st = args$st
#if (st != "yes" || st != "no"){
#print ("Error please provide yes or no for -st")
#quit()
#}

inpsamples = args$inpsamples

rt = list.files(args$inpsamples)
smvec = strsplit(args$samplenames,",")[[1]]

prefix = args$prefix
remprefix = paste(prefix,"_",sep="")
tvec = gsub(remprefix,"",smvec)
flankingbases = args$flankingbases
ft = args$ft

outdir = args$outdir
dir.create(outdir, showWarnings = FALSE)

for (i in (1:length(smvec))){
print (smvec[i])
if (st == "yes"){
scn = rt[grep(paste("neg_",smvec[i],"_scores.txt",sep=""),rt)]
wtn = rt[grep(paste("neg_",smvec[i],"_weight.txt",sep=""),rt)]

scp = rt[grep(paste("pos_",smvec[i],"_scores.txt",sep=""),rt)]
wtp = rt[grep(paste("pos_",smvec[i],"_weight.txt",sep=""),rt)]

print (scn)
datan <- read.delim(file = paste(inpsamples,scn,sep="/"), sep = "\t", h = F)
datan = datan[((datan$V1)>=(-flankingbases)),]
datan = datan[((datan$V1)<=(ft+flankingbases)),]
datan$V1 <- round(datan$V1,0)
datan$V2 <- round(datan$V2,2)
datan = datan[!duplicated(datan), ]
datan = data.frame(datan %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
datan$V2 = na.locf(datan$V2)
datan = datan[,c("V1","V2")]
datan$V1 = (datan$V1-1000)*-1

print (scp)
data <- read.delim(file = paste(inpsamples,scp,sep="/"), sep = "\t", h = F)
data = data[((data$V1)>=(-flankingbases)),]
data = data[((data$V1)<=(ft+flankingbases)),]
data$V1 <- round(data$V1,0)
data$V2 <- round(data$V2,2)
data = data[!duplicated(data), ]
data = data.frame(data %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
data$V2 = na.locf(data$V2)
data = data[,c("V1","V2")]
data = rbind(datan, data)
datas_agg <- aggregate(V2 ~ V1 , data, median)
datas_agg$sample <- rep(tvec[i], nrow(datas_agg))
colnames(datas_agg) <- c("Distance", "likelihood", "sample")

print (wtn)
datawn <- read.delim(file = paste(inpsamples,wtn,sep="/"), sep = "\t", h = F)
datawn = datawn[((datawn$V1)>=(-flankingbases)),]
datawn = datawn[((datawn$V1)<=(ft+flankingbases)),]
datawn$V1 <- round(datawn$V1,0)
datawn = datawn[!duplicated(datawn), ]
datawn = data.frame(datawn %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
datawn$V2 = na.locf(datawn$V2)
datawn = datawn[,c("V1","V2")]
datawn$V1 = (datawn$V1-1000)*-1

print (wtp)
dataw <- read.delim(file = paste(inpsamples,wtp,sep="/"), sep = "\t", h = F)
dataw = dataw[((dataw$V1)>=(-flankingbases)),]
dataw = dataw[((dataw$V1)<=(ft+flankingbases)),]
dataw$V1 <- round(dataw$V1,0)
dataw = dataw[!duplicated(dataw), ]
dataw = data.frame(dataw %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
dataw$V2 = na.locf(dataw$V2)
dataw = dataw[,c("V1","V2")]
dataw = rbind(datawn, dataw)
dataw_agg <- aggregate(V2 ~ V1 , dataw, sum)
dataw_agg$sample <- rep(tvec[i], nrow(dataw_agg))
colnames(dataw_agg) <- c("Distance", "likelihood", "sample")

} else if (st == "no") {

sc = rt[grep(paste(smvec[i],"_scores.txt",sep=""),rt)]
wt = rt[grep(paste(smvec[i],"_weight.txt",sep=""),rt)]

print (sc)
data <- read.delim(file = paste(inpsamples,sc,sep="/"), sep = "\t", h = F)
data = data[((data$V1)>=(-flankingbases)),]
data = data[((data$V1)<=(ft+flankingbases)),]
data$V1 <- round(data$V1,0)
data$V2 <- round(data$V2,2)
data = data[!duplicated(data), ]
data = data.frame(data %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
data$V2 = na.locf(data$V2)
data = data[,c("V1","V2")]
datas_agg <- aggregate(V2 ~ V1 , data, median)
datas_agg$sample <- rep(tvec[i], nrow(datas_agg))
colnames(datas_agg) <- c("Distance", "likelihood", "sample")

print (wt)
dataw <- read.delim(file = paste(inpsamples,wt,sep="/"), sep = "\t", h = F)
dataw = dataw[((dataw$V1)>=(-flankingbases)),]
dataw = dataw[((dataw$V1)<=(ft+flankingbases)),]
dataw$V1 <- round(dataw$V1,0)
dataw = dataw[!duplicated(dataw), ]
dataw = data.frame(dataw %>% complete(nesting(V3), V1 = seq(min(V1), max(V1), 1L)))
dataw$V2 = na.locf(dataw$V2)
dataw = dataw[,c("V1","V2")]
dataw_agg <- aggregate(V2 ~ V1 , dataw, sum)
dataw_agg$sample <- rep(tvec[i], nrow(dataw_agg))
colnames(dataw_agg) <- c("Distance", "likelihood", "sample")

}

data = merge(datas_agg,dataw_agg,by=c("Distance"),all=TRUE)
data[,c(2,4)][is.na(data[,c(2,4)])] = 0
data[,2] = as.numeric(data[,2])
data[,4] = as.numeric(data[,4])
data$wa = data[,2]*data[,4]
data_agg = data[,c(1,6,3)]
colnames(data_agg) = c("Distance","likelihood","sample")

if (i == 1){
datafinal = data_agg
avgdata = datas_agg
countdata = dataw_agg
}
else {
datafinal = rbind(datafinal, data_agg)
avgdata = rbind(avgdata,datas_agg)
countdata = rbind(countdata,dataw_agg)
}
}
finallist = list()
finallist[[1]] = datafinal
finallist[[2]] = avgdata
finallist[[3]] = countdata
save(file = paste(paste(outdir,prefix,sep="/"),"_MetaORF_scores.RData",sep=""), finallist)

print (completed)
