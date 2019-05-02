
# Author:      Mohamood Adhil
# Date:        10th July 2017 
# help:        Rscript <R-script> -h (or) --help

#######
library(argparse)
#######
parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to perform replication-transcription collision analysis using replication origin and transcript directionality')
parser$add_argument("-sc", dest="sc", help="bed file containing gene information with strand information, third column is strand information", default="na", required = TRUE)
parser$add_argument("-rep", dest="rep", help="Replication Origin bed file, if not provided then -arep file should be given", default="na")
parser$add_argument("-arep", dest="arep", help="Annotated Replication Origin bed file, if not provided then -rep should be provided", default="na")
parser$add_argument("-bc", dest="bc", help="comma seperated bed sample file path", required = TRUE)
parser$add_argument("-sn", dest="sn", help="sample name comma seperated for collisiontype", required = TRUE)
parser$add_argument("-pr", dest="pr", required=TRUE, help="Prefix collisiontype for output file (Mandatory)")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

###################################################################

library(bedr)
library(ggplot2)

annotfile = args$sc
repfile = args$rep
arepfile = args$arep
fn = args$bc 
sbdg = strsplit(fn,",")[[1]]
sn = args$sn
stags = strsplit(sn,",")[[1]]
pref = args$pr
outdir = args$dir

if ((repfile == "na")&(arepfile == "na")){
 print("Please provide atleast one file -sc or -dc")
 q()
}  

strandfunc = function(sc,findata){
  sc$gid = paste(sc$gchr,":",sc$gstart,"-",sc$gend,sep="")
  scv = sc$gid
  scv.sort = bedr.sort.region(scv)
  findata$pid = paste(findata$chr,":",findata$start,"-",findata$end,sep="")
  dtv = findata$pid
  dtv.sort = bedr.sort.region(dtv)
  dtv.int <- bedr(input = list(a = dtv.sort, b = scv.sort), engine="/usr/bin/bedtools", method="intersect",params ="-wo -sorted")
  #dtv.int = dtv.int[-which(dtv.int$V5 == "-1"),]
  dtv.int$gid = paste(dtv.int$V4,":",dtv.int$V5,"-",dtv.int$V6,sep="")
  dfin = dtv.int[,c(1,6,5)]
  colnames(dfin) = c("pid","gid","Overlap")
  findata1 = merge(findata,dfin,by=c("pid"))
  findata1 = merge(findata1,sc,by=c("gid"))
  findata1$origid = as.vector(findata1$origid)
  findata1$firing = as.vector(findata1$firing)
  return (findata1)
}

sc = read.delim(file = annotfile,sep="\t",h=F)
colnames(sc) = c("gchr","gstart","gend","gstrand","gname","gsymbol","type")

dlength = c(250,500,1000,2000,5000,10000)
if (repfile != "na"){
reporig = read.delim(file = repfile,sep="\t",h=F,skip=1)
reporig = reporig[,c(1:7)]
colnames(reporig) = c("chromosome","start","end","strand","type","originid","firing")
reporig$center = round(reporig$start + (reporig$end - reporig$start)/2, 0)
reqcat = c("Early/middle","Early")
reporig=reporig[which(reporig$firing %in% reqcat),]
for (i in 1:length(dlength)){
tb = data.frame(chr = reporig$chromosome, start = reporig$center-dlength[i], end = reporig$center, origid = reporig$originid, firing = reporig$firing, typeo="l")
tb2 = data.frame(chr = reporig$chromosome, start = reporig$center+1, end = reporig$center+dlength[i], origid = reporig$originid, firing = reporig$firing, typeo="r")
tbx = rbind(tb,tb2)
tbx$start[which(tbx$start <= 0)] = 0
tbnp = strandfunc(sc,tbx)
tbnp$collision = "crt"
tbnp$typeo = as.vector(tbnp$typeo)
tbnp$gstrand = as.vector(tbnp$gstrand)
tbnp[which(tbnp$typeo == "l" & tbnp$gstrand == "1"),"collision"] = paste("headon",dlength[i],sep="")
tbnp[which(tbnp$typeo == "l" & tbnp$gstrand == "-1"),"collision"] = paste("codirectional",dlength[i],sep="")
tbnp[which(tbnp$typeo == "r" & tbnp$gstrand == "1"),"collision"] = paste("codirectional",dlength[i],sep="")
tbnp[which(tbnp$typeo == "r" & tbnp$gstrand == "-1"),"collision"] = paste("headon",dlength[i],sep="")
tbnp[which(tbnp$typeo == "l" & tbnp$gstrand == "1"),"collisiontype"] = "headon"
tbnp[which(tbnp$typeo == "l" & tbnp$gstrand == "-1"),"collisiontype"] = "codirectional"
tbnp[which(tbnp$typeo == "r" & tbnp$gstrand == "1"),"collisiontype"] = "codirectional"
tbnp[which(tbnp$typeo == "r" & tbnp$gstrand == "-1"),"collisiontype"] = "headon"
tbnp$DistanceGroup = dlength[i]
dups = as.vector(tbnp$gname[duplicated(tbnp$gname)])
tbnp = tbnp[which(!tbnp$gname %in% dups),]
tbnp = tbnp[,c(10,11,12,14,17,18,19)]
if (i == 1){
alldata = tbnp
}else{
alldata = rbind(alldata,tbnp)
}
}
fname = paste(outdir,"allgenes_reptranscollision.txt",sep="/")
write.table(file = fname,alldata,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
alldata$Count = 1 
disfreq = aggregate(Count ~ collision+collisiontype+DistanceGroup, alldata, sum)
disfreq$DistanceGroup = factor(disfreq$DistanceGroup,levels=dlength)
allsum = aggregate(Count ~ DistanceGroup, alldata, sum)
colnames(allsum)[2] = c("TotalGenes")
disfreq = merge(disfreq, allsum, by=c("DistanceGroup"))
disfreq$perc = round(disfreq$Count/disfreq$TotalGenes,2)*100
p <- ggplot(data=disfreq, aes(x=DistanceGroup, y=perc, fill=collisiontype)) + geom_bar(stat="identity", position=position_dodge()) + xlab("Replication Origin to Transcript Collision Distance") + ylab("Percentage of Gene numbers") 
fname = paste(outdir,"ReplicationTranscriptionCollision.jpeg",sep="/")
ggsave(p, file = fname, height = 6, width = 8, dpi = 300)
fname = paste(outdir,"ReplicationTranscriptionCollision_aggregate.txt",sep="/")
write.table(file = fname,disfreq,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}
if(arepfile != "na"){
alldata = read.delim(file = arepfile,sep="\t",h=F,skip=1)
}

for (i in 1:length(sbdg)){
  findata = read.delim(file = sbdg[i],sep= "\t",h=F)
  colnames(findata)=c("chr","start","end","width","score","summit","peakid")
  tbx = strandfunc(alldata,findata)
  tbx$Overlap = as.numeric(tbx$Overlap)
  tbx$gwidth = tbx$gend - tbx$gstart
  tbx$Overlap = as.numeric(tbx$Overlap)
  tbx = unique(tbx)
  tbx$sample = stags[i]
  atbx = aggregate(Overlap ~ gname+gwidth+DistanceGroup+collisiontype+collision, tbx, sum)
  gsum = aggregate(gwidth ~ DistanceGroup+collisiontype+collision, atbx, sum)
  osum = aggregate(Overlap ~ DistanceGroup+collisiontype+collision, atbx, sum)
  fdat = merge(osum, gsum, by=c("collision"))
  fdat$perc = (fdat$Overlap/fdat$gwidth)*100
  fdat$sample = stags[i]
  fdat = fdat[,c(2,3,4,7,8)]
  colnames(fdat) = c("Distance", "Collision", "Overlap", "GeneBase", "Percentage")
  fdat$Distance=factor(fdat$Distance,levels=dlength)
  yname = paste(stags[i],"Percentage of peak base coverage", sep=" ")
  p <- ggplot(data=fdat, aes(x=Distance, y=Percentage, fill=Collision)) + geom_bar(stat="identity", position=position_dodge()) + xlab("Replication Origin to Transcript Collision Distance") + ylab(yname)
  fname = paste(outdir,"/",stags[i],"_RepTrans_Collision_Percentage.jpeg",sep="")
  ggsave(p, file = fname, height = 6, width = 8, dpi = 300)
  fdat$sample = stags[i]
  if (i == 1){
  alltbx = tbx
  allfdat = fdat	
  }else{
  alltbx = rbind(alltbx, tbx)
  allfdat = rbind(allfdat, fdat)
  }
}
alltbx = alltbx[,c(3:7,10,12:14,16,17,20)]
colnames(alltbx)[7:9] = c("genestart","geneend","geneID")
fname = paste(outdir,"/",pref,"_RepTrans_Collision_allgenes.txt",sep="")
write.table(file = fname,alltbx,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
fname = paste(outdir,"/",pref,"_RepTrans_Collision_aggregate.txt",sep="")
write.table(file = fname,allfdat,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)



