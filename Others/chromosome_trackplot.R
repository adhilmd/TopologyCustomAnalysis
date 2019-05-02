# Author:      Mohamood Adhil
# Date:        15th July 2017 
# help:        Rscript <R-script> -h (or) --help

#######
library(argparse)
#######
parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to plot the user defined regions using the bed and bed graph files')
parser$add_argument("-chr", dest="chr", help="chromosome number to plot (Mandatory)", required = TRUE)
parser$add_argument("-st", dest="st", type="integer", help="start site (Mandatory)", required = TRUE)
parser$add_argument("-en", dest="en", type="integer", help="end site (Mandatory)", required = TRUE)
parser$add_argument("-bdg", dest="bdg", help="Comma seperated bedgraph files, only four allowed (Mandatory)", required = TRUE)
parser$add_argument("-bd", dest="bd", help="Comma seperated bed files, only four allowed and also should match the number of bedgraph files (Mandatory)", required = TRUE)
parser$add_argument("-sc", dest="sc", required=TRUE, help="Protein coding bed annotation file (Mandatory)")
parser$add_argument("-tg", dest="tg", required=TRUE, help="Comma seperated tags, should match the number of bed and bedgraph files (Mandatory)")
parser$add_argument("-pr", dest="pr", required=TRUE, help="Prefix for output file (Mandatory)")
parser$add_argument("-ma", dest="ma", help="Main Title for the plot", default="")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()


###################################################################

library(Gviz)
library(GenomicRanges)

chromnum = args$chr
startsite = args$st
endsite = args$en

bg = args$bdg
tg = args$tg
bedf = args$bd

bgl = strsplit(bg,",")[[1]] 
tgl = strsplit(tg,",")[[1]]
bdl = strsplit(bedf,",")[[1]]

if (length(bgl) != length(bdl)){
print ("Error: Provide the same number of files in -bdg and -bd")
q()
}

nsample = length(bgl)

pcoding = args$sc

prefix = args$pr

maintitle = args$ma

outdir = args$dir

allcolors = c("orange","violetred","midnightblue","red")

chrc = function(dat){
  chrm = c("chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9","plasmid")
  chrr = c("X","XI","XII", "XIII", "XIV", "XV", "XVI", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX","plasmid")
  dat$seqnames = as.vector(dat$seqnames)
  for (i in 1:length(chrm)){
    dat[which(dat$seqnames == chrm[i]),1] = chrr[i] 
  }
  return(dat)
}

options(ucscChromosomeNames=FALSE)

pc = read.delim(file = pcoding,sep="\t",h=F,skip=1)
pc = pc[,c(1:7)]
pc$V8 = pc$V3-pc$V2
pc$V9 = "-"
pc$V10 = "-"
pc = pc[,c(1,2,3,8,4,7,5,6)]
colnames(pc)[1] = "seqnames"
pc = chrc(pc)
colnames(pc) = c("chromosome","start","end","width","strand","feature","gene","symbol")
pc$feature = as.vector(pc$feature)
pc$symbol = as.vector(pc$symbol)
pc$gene = as.vector(pc$gene)
pc$symbol[which(pc$symbol == "NULL")] <- pc$gene[which(pc$symbol == "NULL")]
grtr <- GeneRegionTrack(pc, name="Genes", strand = pc$strand, showId = TRUE, cex.axis=1)
AT=GenomeAxisTrack()

dtrack = list()
nptrack = list()
for (i in 1:nsample){
data1 = read.delim(file = bgl[[i]][1],sep="\t",h=F)
data1$strand="*"
colnames(data1) = c("seqnames","start","end","score","strand")
data1 = data1[,c("seqnames","start","end","strand","score")]
data1$strand = "*"
data1$score = as.numeric(as.vector(data1$score))
data1 = chrc(data1)
gfdata1 <- makeGRangesFromDataFrame(data1, keep.extra.columns=TRUE)
data1.track = DataTrack(gfdata1,strand="*",genome="Sc03",col.histogram=allcolors[i],fill.histogram=allcolors[3], name=tgl[[i]], col.sampleNames="black",col.title="black",col.axis="black",cex.axis=0.5,fontfamily.legend="arial")
dtrack[[i]] = data1.track

allbed = read.delim(file = bdl[[i]][1],sep="\t",h=F)
allbed = allbed[,c(1:3)]
colnames(allbed) = c("seqnames","start","end")
allbed$strand = "*"
allbed = chrc(allbed)
gfbed <- makeGRangesFromDataFrame(allbed, keep.extra.columns=TRUE)
atr.track <- AnnotationTrack(gfbed, name=tgl[[i]],cex.axis=0.5)
nptrack[[i]] = atr.track
}

if (nsample == 1){
  ofile = paste(paste(outdir,prefix,sep="/"),"genes_track.jpeg",sep="_")
  print (ofile)
  jpeg(ofile,width=12,height=8,res=1000,units="in",family="arial")
  plotTracks(c(nptrack[[1]],dtrack[[1]],grtr,AT),chromosome=chromnum,from=startsite,to=endsite,transcriptAnnotation='symbol',window='auto',type="histogram", cex.title=1.5, cex.main = 1.5, main=maintitle, col.main="black")
  dev.off()
}

if (nsample == 2){
  ofile = paste(paste(outdir,prefix,sep="/"),"genes_track.jpeg",sep="_")
  print (ofile)
  jpeg(ofile,width=12,height=8,res=500,units="in",family="arial")
  plotTracks(c(nptrack[[1]],dtrack[[1]],nptrack[[2]],dtrack[[2]],grtr,AT),Neg="darkred",Pos="darkgreen",Stable="blue",chromosome=chromnum,from=startsite,to=endsite,transcriptAnnotation='symbol',window='auto',type="histogram", cex.title=1.5, cex.main = 1.5, main=maintitle, col.main="black")
  dev.off()
}

if (nsample == 3){
  ofile = paste(paste(outdir,prefix,sep="/"),"genes_track.jpeg",sep="_")
  print (ofile)
  jpeg(ofile,width=12,height=8,res=500,units="in",family="arial")
  plotTracks(c(nptrack[[1]],dtrack[[1]],nptrack[[2]],dtrack[[2]],nptrack[[3]],dtrack[[3]],grtr,AT),chromosome=chromnum,from=startsite,to=endsite,transcriptAnnotation='symbol',window='auto',type="histogram", cex.title=1.5, cex.main = 1.5, main=maintitle, col.main="black")
  dev.off()
}

if (nsample == 4){
  ofile = paste(paste(outdir,prefix,sep="/"),"genes_track.jpeg",sep="_")
  print (ofile)
  jpeg(ofile,width=12,height=8,res=500,units="in",family="arial")
  plotTracks(c(nptrack[[1]],dtrack[[1]],nptrack[[2]],dtrack[[2]],nptrack[[3]],dtrack[[3]],nptrack[[4]],dtrack[[4]],grtr,AT),chromosome=chromnum,from=startsite,to=endsite,transcriptAnnotation='symbol',window='auto',type="histogram", cex.title=1.5, cex.main=1.5, main=maintitle, col.main="black")
  dev.off()
}
