
# Author:      Mohamood Adhil
# Date:        24th May 2017 
# help:        Rscript <R-script> -h (or) --help

#######
library(argparse)
#######
parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is to perform intergenic spaces analysis with respect to gene direction of upstream and downstream')
parser$add_argument("-sc", dest="sc", help="bed file containing gene information with strand information, third column is strand information, if not provided then -dc should be provided", default="na")
parser$add_argument("-dc", dest="dc", help="Directionality annotated file, if not provided then -sc should be provided", default="na")
parser$add_argument("-bc", dest="bc", help="comma seperated bed sample file path (Mandatory)", required = TRUE)
parser$add_argument("-sn", dest="sn", help="sample name comma seperated for tag (Mandatory)", required = TRUE)
parser$add_argument("-pr", dest="pr", required=TRUE, help="Prefix tag for output file (Mandatory)")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

###################################################################

library(bedr)
library(ggplot2)

annotfile = args$sc
directionalityfile = args$dc
fn = args$bc 
fpath = strsplit(fn,",")[[1]]
sn = args$sn
sname = strsplit(sn,",")[[1]]
pref = args$pr
outdir = args$dir

if ((annotfile == "na")&(directionalityfile == "na")){
 print("Please provide atleast one file -sc or -dc")
 q()
}  

slit <- function(sc,chr){
  alist = list()
  j = 0
  for (i in chr){
    j = j+1
    te = sc[which(sc$V1 == i),]
    te = te[order(te$V2,te$V3),]
    alist[[j]] = te
  }
  adf = do.call("rbind",alist)
  return(adf)
}

if (annotfile != "na"){
sc = read.delim(file = annotfile, sep ="\t", h = F, skip = 1)
sc = sc[,c(1:7)]
sc$V1 = as.vector(sc$V1)
chr = paste("chr",c(1:16),sep="")
adf = slit(sc,chr)
bdf = adf[-1,]
colnames(bdf) = paste("V",c(8:14),sep="")
adf = adf[-dim(adf)[1],]
cdf = cbind(adf,bdf)
cdf$V4 = as.vector(cdf$V4)
cdf$V11 = as.vector(cdf$V11)
tagsdirc = c()
for (i in 1:dim(cdf)[1]){
  if (cdf[i,4] == "1" && cdf[i,11] == "1"){
  tagsdirc = c(tagsdirc,"Codirectional+")
  }
  if (cdf[i,4] == "-1" && cdf[i,11] == "-1"){
    tagsdirc = c(tagsdirc,"Codirectional-")
  }
  if (cdf[i,4] == "-1" && cdf[i,11] == "1"){
    tagsdirc = c(tagsdirc,"Divergent")
  }
  if (cdf[i,4] == "1" && cdf[i,11] == "-1"){
    tagsdirc = c(tagsdirc,"Convergent")
  }
}
cdf$V15 = tagsdirc
cdf = cdf[-which(cdf$V1 != cdf$V8),]
dat1 = cdf
dat1 = dat1[-which(dat1$V3 > dat1$V9),]
dat1$V16 = dat1$V1
dat1$V17 = dat1$V3
dat1$V18 = dat1$V9
dat1$len1 = dat1$V3-dat1$V2
dat1$len2 = dat1$V10-dat1$V9
dat1 = dat1[,c(16,17,18,15,5,12,19,20)]
datall = dat1
datall = datall[which(datall$V17 != datall$V18),]
datall$id = paste("p",c(1:dim(datall)[1]),sep="")
datall$length = datall$V18-datall$V17
colnames(datall) = c("gchr","fend","sstart","tag","fgeneid","sgeneid","len1","len2","id","distance")
tdis = datall[order(datall$distance),]
t1 = tdis[which(tdis$distance >=1 & tdis$distance <= 250),]
t1$cdis = "<=250"
t2 = tdis[which(tdis$distance >=251 & tdis$distance <= 500),]
t2$cdis = "251-500"
t3 = tdis[which(tdis$distance >501),]
t3$cdis = ">500"
directiondat = rbind(t1,t2,t3)
directiondat$cdis = factor(directiondat$cdis,levels=c("<=250","251-500",">500"))
write.table(file = paste(outdir,"Directionality_proteincoding.txt",sep="/"), directiondat, sep ="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

disfreq = data.frame(table(directiondat$tag, directiondat$cdis))
colnames(disfreq) = c("Tag","DistanceGroup","Count")
disfreq$DistanceGroup = factor(disfreq$DistanceGroup,levels=c("<=250","251-500",">500"))
colnames(disfreq)[1] = c("IntergenicSpaces")
p <- ggplot(data=disfreq, aes(x=DistanceGroup, y=Count, fill=IntergenicSpaces)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Intergenic Distance") + ylab("Number of Gene Pairs") 
fname = paste(outdir,"Genearrangement.jpeg",sep="/")
ggsave(p, file = fname, height = 6, width = 8, dpi = 300)
}

if (directionalityfile != "na"){
directiondat = read.delim(file = directionalityfile, sep ="\t", h = T)
}

mefunc = function(sc,findata){
  sc$gid = paste(sc$gchr,":",sc$fend,"-",sc$sstart,sep="")
  scv = sc$gid
  scv.sort = bedr.sort.region(scv)
  findata$pid = paste(findata$chromosome,":",findata$start,"-",findata$end,sep="")
  dtv = findata$pid
  dtv.sort = bedr.sort.region(dtv)
  dtv.int <- bedr(input = list(a = dtv.sort, b = scv.sort), engine="bedtools", method="intersect",params ="-wo -sorted")
  dtv.int$gid = paste(dtv.int$V4,":",dtv.int$V5,"-",dtv.int$V6,sep="")
  dfin = dtv.int[,c(1,6,5)]
  colnames(dfin) = c("pid","gid","Overlap")
  findata1 = merge(findata,dfin,by=c("pid"))
  findata1 = merge(findata1,sc,by=c("gid"),all=TRUE)
  findata1$fgeneid = as.vector(findata1$fgeneid)
  findata1$tag = as.vector(findata1$tag)
  findata1$Overlap = as.numeric(as.vector(findata1$Overlap))
  findata1$distance = as.numeric(as.vector(findata1$distance))
  return (findata1) 
}

for (f in 1:length(fpath)){
findata = read.delim(file = fpath[f],sep= "\t",h=F)
colnames(findata) = c("chromosome","start","end","length","score","summit","peakid")
findata1 = mefunc(directiondat,findata) 
findata1[is.na(findata1$Overlap),"Overlap"] = 0
findata1 = findata1[,c(14,20,10,21)]
odissum = aggregate(distance ~ tag+cdis, directiondat, sum)
overlapsum = aggregate(Overlap ~ tag+cdis, findata1, sum)
fdat = merge(odissum, overlapsum, by=c("tag","cdis"))
fdat['perc']=(fdat['Overlap']/fdat['distance'])*100
fdat$sample = sname[f]
colnames(fdat)[1] = c("IntergenicSpaces")
yname = paste(sname[f],"Base Percentage", sep=" ")
p <- ggplot(data=fdat, aes(x=cdis, y=perc, fill=IntergenicSpaces)) + geom_bar(stat="identity", position=position_dodge()) + xlab("Intergenic Distance") + ylab(yname)
fname = paste(outdir,"/",sname[f],"_directionalityIntergenicspaces.jpeg",sep="")
ggsave(p, file = fname, height = 6, width = 8, dpi = 300)
if (f == 1){
alldata = fdat
}
else{
alldata = rbind(alldata,fdat)
}
fname = paste(outdir,"/",pref,"_directionalityIntergenicspaces.txt",sep="")
write.table(file = fname, alldata, sep ="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
