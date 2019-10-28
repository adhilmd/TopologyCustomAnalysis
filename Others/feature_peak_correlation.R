# Author:      Mohamood Adhil
# Date:        25th July 2017 
# help:        Rscript <R-script> -h (or) --help

#######
library(argparse)
#######
parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to perform correlation analysis between two user bed files with respect to gene coverage')
parser$add_argument("-sc", dest="sc", help="bed file containing gene information with strand information, third column is strand information (Mandatory)", required = TRUE)
parser$add_argument("-ud", dest="ud", type="integer", help="Upstream and Downstream flanking bases to include", default=0)
parser$add_argument("-bd", dest="bd", help="comma seperated two bed files path for sample1VSsample2 (Mandatory)", required = TRUE)
parser$add_argument("-tg", dest="tg", help="comma seperated tags, should be 2 tags sample1 and sample2 (Mandatory)", required = TRUE)
parser$add_argument("-pr", dest="pr", required=TRUE, help="Prefix for output file (Mandatory)")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

###################################
library(ggplot2)
library(bedr)

options("scipen"=100, "digits"=4)

annotfile = args$sc

fname = args$bd
bfiles = strsplit(fname,",")[[1]]

tgs = args$tg
atags = strsplit(tgs,",")[[1]]

udb = args$ud

pref = args$pr

outdir = args$dir
dir.create(outdir, showWarnings = FALSE)

#bfiles = c("/lustre/home/amohamme/data/Chip/NegSupercoiling/R1308/Output/DataProcessing/bed_graphfiles/proc_R1308_WT_G1_IP-A.bed","/lustre/home/amohamme/data/Chip/NegSupercoiling/R1308/Output/DataProcessing/bed_graphfiles/proc_R1308_WT_S_IP-A.bed")

mefunc = function(sc,findata){
  sc$gid = paste(sc$gchr,":",sc$gsu,"-",sc$ged,sep="")
  scv = sc$gid
  scv.sort = bedr.sort.region(scv)
  findata$pid = paste(findata$chr,":",findata$start,"-",findata$end,sep="")
  dtv = findata$pid
  dtv.sort = bedr.sort.region(dtv)
  dtv.int <- bedr(input = list(a = dtv.sort, b = scv.sort), engine="bedtools", method="intersect",params ="-wo -sorted")
  #dtv.int = dtv.int[-which(dtv.int$V5 == "-1"),]
  dtv.int$gid = paste(dtv.int$V4,":",dtv.int$V5,"-",dtv.int$V6,sep="")
  dfin = dtv.int[,c(1,6,5)]
  colnames(dfin) = c("pid","gid","Overlap")
  findata1 = merge(findata,dfin,by=c("pid"))
  findata1 = merge(findata1,sc,by=c("gid"))
  findata1$geneid = as.vector(findata1$geneid)
  findata1$Type = as.vector(findata1$Type)
  findata1$Overlap = as.numeric(as.vector(findata1$Overlap))
  findata1$gwidth = as.numeric(as.vector(findata1$gwidth))
  findata1$width = as.numeric(as.vector(findata1$width))
  return (findata1) 
}

sc = read.delim(file = annotfile,sep="\t",h=F,skip=1)
sc = sc[,c(1:7)]
colnames(sc) = c("gchr","gstart","gend","gstrand","geneid","genesymbol","genetype")
sc$gwidth = sc$gend-sc$gstart
sc = sc[-which(sc$gwidth == 0),]
sc$gsu = sc$gstart - udb
sc$ged = sc$gend + udb
sc$gsu[which(sc$gsu <= 0)] = 0
sc$gewidth = sc$ged - sc$gsu

mfunc = function(bfile,sc){
  findata = read.delim(file = bfile, sep="\t", h=F)
  colnames(findata) = c("chr","start","end","width","score","summit","ID")
  coild = mefunc(sc,findata)
  coild = coild[,c(3:20)]
  colnames(coild) = c("chr","Peak_start","Peak_end","width","score","summit","id","Overlap","Gene_Chr","Gene_Start","Gene_End","strand","Gene_id","Gene_symbol","Feature","Gene_Length","gsu","ged")
  ind = which(coild$Gene_Length == 0)
  if (length(ind) > 0){
    coild = coild[-ind,]
  }
  coild$alllength = coild$ged-coild$gsu
  reqdata = aggregate(Overlap~Gene_id+alllength, coild, sum)
  colnames(reqdata) = c("geneid","length","udoverlap")
  reqdata$perc = reqdata$udoverlap/reqdata$length
  return (reqdata)
}

sample1 = mfunc(bfiles[1],sc)
sample2 = mfunc(bfiles[2],sc)

mesample = merge(sample1, sample2, by = c("geneid"))
mesample = na.omit(mesample)
mesample$perc.x = (mesample$perc.x)*100
mesample$perc.y = (mesample$perc.y)*100

print("T-TEST")
t.test(mesample$perc.x,mesample$perc.y,paired=TRUE)


corresult = cor.test(mesample$perc.x,mesample$perc.y)
print (corresult)
stat = data.frame(c(corresult$estimate,corresult$p.value,corresult$statistic))
row.names(stat) = c("Correlation","P-value","t-statistic") 
colnames(stat) = "Values"
fname = paste(outdir,"/",pref,"_correlation.txt",sep="")
write.table(file = fname, stat, row.names = TRUE, col.names = FALSE, quote = FALSE)

xtag = atags[1]
ytag = atags[2]

fname = paste(outdir,"/",pref,"_correlation.jpeg",sep="")
p = ggplot(data = mesample, aes(x = perc.x, y = perc.y)) + geom_point(color='red') + geom_smooth(method = "lm", se = FALSE) + xlab(xtag) + ylab(ytag)
ggsave(p, file = fname, height = 6, width = 8, dpi = 300)
