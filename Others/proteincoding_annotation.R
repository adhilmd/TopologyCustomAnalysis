# Author:      Mohamood Adhil
# Date:        18th April 2017 
# help:        Rscript <R-script> -h (or) --help

#######
library(argparse)
#######
parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to plot the user defined regions using the bed and bed graph files')
parser$add_argument("-sc", dest="sc", help="Protein coding annotation bed file (Mandatory)", required = TRUE)
parser$add_argument("-bd", dest="bd", help="bed file (Mandatory)", required = TRUE)
parser$add_argument("-udb", dest="udb", type="integer", help="Upstream and Downstream bases for annotation (Mandatory)", required = TRUE)
parser$add_argument("-pr", dest="pr", required=TRUE, help="Prefix for output file (Mandatory)")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

####################################################################################

library(bedr)

options("scipen"=100, "digits"=4)

bed = args$bd
pref = args$pr
annotfile = args$sc
udb = args$udb

outdir = args$dir
dir.create(outdir, showWarnings = FALSE)

sc = read.delim(file = annotfile,sep="\t",h=F,skip=1)
sc = sc[,c(1:7)]
colnames(sc) = c("gchr","gstart","gend","gstrand","geneid","genesymbol","genetype")
sc$ustart = sc$gstart - udb
sc$dend = sc$gend + udb
sc$ustart[which(sc$ustart <= 0)] = 0
sc$gwidth = sc$gend-sc$gstart
sc = sc[,c(1,8,9,2,3,4,5,6,7,10)]

mefunc = function(sc,findata){
  sc$gid = paste(sc$gchr,":",sc$ustart,"-",sc$dend,sep="")
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
  findata1$Type = as.vector(findata1$Type)
  findata1$width = as.numeric(as.vector(findata1$width))
  return (findata1) 
}

findata = read.delim(file = bed, sep="\t", h=F) 
colnames(findata) = c("chr","start","end","width","score","summit","ids")
datfin = mefunc(sc,findata)

pfin = datfin[which(datfin$gstrand == "1"),]
nfin = datfin[which(datfin$gstrand == "-1"),]

upstream = pfin[which(pfin$start >= pfin$ustart & pfin$start < pfin$gstart & pfin$end > pfin$ustart & pfin$end <= pfin$gstart),]
temp = pfin[which(pfin$start < pfin$ustart & pfin$start < pfin$gstart & pfin$end > pfin$ustart & pfin$end <= pfin$gstart),]
upstream = rbind(upstream,temp)
if (dim(upstream)[1] >= 1){
upstream$Annottype = paste("Upstream",udb,sep="")
}
downstream = pfin[which(pfin$start >= pfin$gend & pfin$start < pfin$dend & pfin$end > pfin$gend & pfin$end <= pfin$dend),]
temp = pfin[which(pfin$start >= pfin$gend & pfin$start < pfin$dend & pfin$end > pfin$gend & pfin$end > pfin$dend),]
downstream = rbind(downstream,temp)
if (dim(downstream)[1] >= 1){
downstream$Annottype = paste("Downstream",udb,sep="")
}

TSS = pfin[which(pfin$start < pfin$gstart & pfin$start < pfin$gend & pfin$end >= pfin$gstart & pfin$end <= pfin$gend),]
TSS$Annottype = "TSS"
TTS = pfin[which(pfin$start >= pfin$gstart & pfin$start <= pfin$gend & pfin$end > pfin$gstart & pfin$end > pfin$gend),]
TTS$Annottype = "TTS"

panot = rbind(upstream,downstream,TSS,TTS)

downstream = nfin[which(nfin$start >= nfin$ustart & nfin$start < nfin$gstart & nfin$end > nfin$ustart & nfin$end <= nfin$gstart),]
temp = nfin[which(nfin$start < nfin$ustart & nfin$start < nfin$gstart & nfin$end > nfin$ustart & nfin$end <= nfin$gstart),]
downstream = rbind(downstream,temp)
if (dim(downstream)[1] >= 1){
downstream$Annottype = paste("Downstream",udb,sep="")
}
upstream = nfin[which(nfin$start >= nfin$gend & nfin$start < nfin$dend & nfin$end > nfin$gend & nfin$end <= nfin$dend),]
temp = nfin[which(nfin$start >= nfin$gend & nfin$start < nfin$dend & nfin$end > nfin$gend & nfin$end > nfin$dend),]
upstream = rbind(upstream,temp)
if (dim(upstream)[1] >= 1){
upstream$Annottype = paste("Upstream",udb,sep="")
}

TTS = nfin[which(nfin$start < nfin$gstart & nfin$start < nfin$gend & nfin$end >= nfin$gstart & nfin$end <= nfin$gend),]
TTS$Annottype = "TTS"
TSS = nfin[which(nfin$start >= nfin$gstart & nfin$start <= nfin$gend & nfin$end > nfin$gstart & nfin$end > nfin$gend),]
TSS$Annottype = "TSS"

nanot = rbind(upstream,downstream,TSS,TTS)

orf = datfin[which(datfin$start >= datfin$gstart & datfin$start < datfin$gend & datfin$end > datfin$gstart & datfin$end <= datfin$gend),]
orf$Annottype = "ORF"
includef = datfin[which(datfin$start < datfin$gstart & datfin$start < datfin$gend & datfin$end > datfin$gstart & datfin$end > datfin$gend),]
includef$Annottype = "IncludeFeature"

finannot = rbind(panot,nanot,orf,includef)
fname = paste(outdir,"/",pref,"_annotation.bed",sep="")
write.table(file = fname, finannot[,c(3:21)], sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

nfeat = data.frame(table(finannot$Annottype))

slices <- c(nfeat$Freq) 
lbls <- as.vector(nfeat$Var1)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls,"-",slices)
lbls <- paste(lbls,"-",pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
fname = paste(outdir,"/",pref,"_annotation.jpeg",sep="")
jpeg(fname, height=2000, width=2400, res=300)
pie(slices,labels = lbls, col=rainbow(length(lbls)), main="", cex=0.8)
dev.off()

#
