setwd("/Users/xung0001/Desktop")
## loading PCC of codng-nonnon group
non <- read.table("cor_noncodingwithrootinfo.txt", sep=" ", header=T)
## loading sense Counts
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
slength <- s$end-s$start
s <- s[!s$chrom=="chloroplast",]
s <- s[!s$chrom=="mitochondria",]
s<- s[,1:8]
## add promoter distance 
non <- non[order(non$cor,decreasing = T),]
non$promoterdis<- c(NA)
length <- c()
for (i in row.names(non))
{print(i)
  if(s[s$name==non[i,2],]$strand=="-")
  {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=s[which(s$name==non[i,2]),]$end-s[which(s$name==non[i,1]),]$start))}
  if(s[s$name==non[i,2],]$strand=="+") {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=s[which(s$name==non[i,1]),]$end-s[which(s$name==non[i,2]),]$start))
  }
  }
## PAT not alway works antisense$end-sense$start
non$promoterdis <- length$length

## add gene length
non$length1 <- c(NA)
non$length2 <- c(NA)
length1 <- c()
length2 <- c()
for (i in row.names(non))
{print(i)
  length1<-rbind(length1,s[which(s$name==non[i,1]),]$end-s[which(s$name==non[i,1]),]$start)
  length2<-rbind(length2,s[which(s$name==non[i,2]),]$end-s[which(s$name==non[i,2]),]$start)
}
non$length1 <- length1
non$length2 <- length2
## adding overlap length
non <- non[order(non$cor,decreasing = T),]
non$overlap <- c(NA)
overlap<- c()
for (i in row.names(non))
{ print(i)
  a1<- s[s$name==non[i,2],]$start 
  a2 <-  s[s$name==non[i,2],]$end 
  b1 <-  s[s$name==non[i,1],]$start
  b2 <- s[s$name==non[i,1],]$end 
  if (a2<=b2 & a1>=b1 & a2>b1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="A",overlap=a2-a1))}
  if (a2>=b2 & a1<=b1 & b1<a2) {overlap<- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="B",overlap=b2-b1))}
  if (a2>b2 & a1<b2 & b1<a1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="C",overlap=b2-a1))}
  if (a2<b2 & a1<b1 & b1<a2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="D",overlap=a2-b1))}
}

overlap <- unique(overlap)
row.names(overlap) <- overlap$id1 
## check where are repeats
non$overlap <- overlap$overlap
overlap[overlap$id1=="AT1G05560",]
overlap[overlap$id1=="AT1G05560",]$overlapc <- "E"
overlap <- overlap[-3,]
overlap[overlap$id1=="AT3G21780",]
overlap <- overlap[-9,]
overlap[overlap$id1=="AT3G21780",]$overlapc <- "E"
## AT4G26510 AT4G26488 divergent pairs AT4G27860 at4g27852 
non$overlapc <- c(NA)
non$overlapc <- overlap$overlapc
non$overlap <- overlap$overlap
## check which ones have strong overlap, we cluter as class E and in this class, overlaplength and promoter distance are very close
non[112,]$overlapc <- "E"
non[non$id1=="AT2G44260",]$overlapc <- "E"
non[non$id1=="AT1G06000",]$overlapc <- "E"
non[non$id1=="AT5G38210",]$overlapc <- "E"
non[non$id1=="AT5G38212",]$overlapc <- "E"
non[non$id1=="AT4G38250",]$overlapc <- "E"
non[non$id1=="AT5G01210",]$overlapc <- "E"
non[non$id1=="AT4G40080",]$overlapc <- "E"
non[non$id1=="AT5G01540",]$overlapc <- "E"
non[non$id1=="AT5G01542",]$overlapc <- "E"
non[non$id1=="AT4G39840",]$overlapc <- "E"
non[non$id1=="AT1G76880",]$overlapc <- "E"
non[non$id1=="AT1G78270",]$overlapc <- "E"
non[non$id1=="AT5G10140",]$overlapc <- "E"
non[non$id1=="AT2G20500",]$overlapc <- "E"
non[non$id1=="AT1G23050",]$overlapc <- "E"
non[non$id1=="AT2G45310",]$overlapc <- "E"
non[non$id1=="AT1G17230",]$overlapc <- "E"
non[non$id1=="AT3G21780",]$overlapc <- "E"
non[non$id1=="AT4G38550",]$overlapc <- "E"
non <- non[,-7]
# load the smedian value for all genes
# load raw counts data
#load counts 
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
s <- s[!s$chrom=="chloroplast",]
slength <- s$end-s$start
row.names(s) <- s$name
s <- s[,-c(1:8)]
# load sample accessions
samples <- as.character(read.table("bams.txt", sep=" ", quote="", header=F, stringsAsFactors=F)[1,])
samples <- gsub("\\.bam","",gsub(".*\\/","",samples))
colnames(s) <- samples
s <- s[,-69]
s <- as.matrix(s)

##normalize data by scales function scale()
s_scale<- t(scale(t(s)))
s_scale <- as.matrix(s_scale)
# normalize gene expression
rpkm<- function(counts,genelength){
  rate <- counts/genelength
  t(t(rate)/colSums(counts))*10^9
}
sr<- rpkm(s,slength)
sr <- na.omit(sr)
s.median_non <- c()
for (i in row.names(non))
{print(i)
  s.median_non <-rbind(s.median_non,data.frame(id1=non[i,1],id2=non[i,2], 
                                               s.median1=median(sr[which(row.names(sr)==non[i,1]),]),
                                               s.median2=median(sr[which(row.names(sr)==non[i,2]),])))
}
s.median_non<- na.omit( s.median_non)
non$s.median <- s.median_non$s.median1
non$as.median <- s.median_non$s.median2
## save the file
write.table(non, file="cor_noncoding_updata.txt", sep="\t", quote=F, row.names=F, col.names=T)
nonc <- read.table("cor_noncoding_updata.txt",header=T)
nonc$overlapc <- as.factor(nonc$overlapc)
summary(nonc$overlapc)
