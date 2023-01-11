##generate correlation value for coding/noncoding groups across 144 Arabidopsis accessions :https://doi.org/10.1038/nature11968
getwd()
setwd("/Users/xung0001/Desktop")
# load raw counts data
#load counts 
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
as <- read.table("counts.antisense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
s <- s[!s$chrom=="chloroplast",]
as <- as[!as$chrom=="chloroplast",]
slength <- s$end-s$start
row.names(s) <- s$name
ann <- s[,1:8]
s <- s[,-c(1:8)]
as <- as[,-c(1:8)]
# load sample accessions
samples <- as.character(read.table("bams.txt", sep=" ", quote="", header=F, stringsAsFactors=F)[1,])
samples <- gsub("\\.bam","",gsub(".*\\/","",samples))
colnames(s) <- samples
#need to remove accession [,69],three ways to perform quality control: A: check the gene counts of seven houskeeping genes: TUB2, TUB8, ACTIN2, UBC1,UBQ10,GAPC, EF-1A in all accessions after normalization in different ways; B: the distribution of raw counts in different accessions; C: raw counts less than 3 will be considered as false positive
s <- s[,-69]
as <- as[,-69]
s <- as.matrix(s)
as <- as.matrix(as)

##normalize data by scales function scale()
s_scale<- t(scale(t(s)))
s_scale <- as.matrix(s_scale)
as_scale<- t(scale(t(as)))
as_scale <- as.matrix(as_scale)
# normalize gene expression
rpkm<- function(counts,genelength){
  rate <- counts/genelength
  t(t(rate)/colSums(counts))*10^9
}
sr<- rpkm(s,slength)
sr <- na.omit(sr)
## load sense-antisense pairs
co_noncoding <- read.table("coding-noncoding.txt", sep=",", header = T)
# select s.median of coding-coding genes greater than 3
s.median <- c()
for (i in row.names(co_noncoding))
{print(i)
  s.median <-rbind(s.median,data.frame(id1=co_noncoding[i,1],id2=co_noncoding[i,2], s.median1=median(s[which(row.names(s)==co_noncoding[i,1]),]),
                                       s.median2=median(s[which(row.names(s)==co_noncoding[i,2]),])))  
}

s.median <- na.omit(s.median)
s.median <-  s.median[which(s.median$s.median1>=3),]
s.median <-  s.median[which(s.median$s.median2>=3),]

### correlation analysis for different groups coding_noncoding
ctr <- data.frame(NULL)
for (i in row.names(s.median))
{print(i)
  ct <- cor.test(s_scale[s.median[i,1],],s_scale[s.median[i,2],])
  mean1 <- mean(sr[which(row.names(sr)==s.median[i,1]),])
  mean2 <- mean(sr[which(row.names(sr)==s.median[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=s.median[i,1], id2=s.median[i,2], cor=ct$estimate,p.value=ct$p.value,
                   stringsAsFactors=F, check.names=F ,mean1=mean1,mean2=mean2, ratio=mean2/mean1
    ))}
hist(ctr$p.value,fre=F)
# notice p.value has been adjusted automatically by cor.test
ctr <- ctr[order(ctr$p.value,decreasing = F),]
# set cut_off to keep data 0.05/number of test=166
dim(ctr)[1]
ctr <- ctr[ctr$p.value<=0.05/dim(ctr)[1],]
ctr <- ctr[order(ctr$cor,decreasing = T),] 
## extract overalpping information for coding-coding 
s1 <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
s1 <- s1[!s1$chrom=="chloroplast",]
s1 <- s1[,-69]
row.names(s1) <- s1$name
noncoding <- ctr
noncoding$length <-

sno1 <- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1<=b1 & 
      a2>=b2
      &
  b1-b2>a1-a2)
  { sno1 <-rbind(sno1,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=(a1-a2)/(b1-b2),overlap=a1-a2)
  )}
}

sno2 <- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if (a2<=b2 & a1>=b1 & a1-a2 >b1-b2)
  { sno2 <-rbind(sno2,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=(b1-b2)/(a1-a2),overlap=b1-b2)
  )}
}
sno3 <- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1<b1 & 
      a2<b2)
    if (a1-a2>b1-b2)
  { sno3 <-rbind(sno3,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(a1-b2)/(b1-b2+a1-a2), overlap=a1-b2)
  )}
}

sno4<- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1>b1 & 
      a2>b2 & a1-a2<b1-b2)
  { sno4 <-rbind(sno4,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(b1-a2)/(b1-b2+a1-a2), overlap=b1-a2)
  )}
}
sno5<- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1<b1 &
      a2<b2 & a1-a2<=b1-b2)
  { sno5 <-rbind(sno5,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(a1-b2)/(b1-b2+a1-a2), overlap=a1-b2)
  )}
}
sno6<- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1>b1 & 
      a2>b2 & a1-a2>=b1-b2)
  { sno6 <-rbind(sno6,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(b1-a2)/(b1-b2+a1-a2), overlap=b1-a2)
  )}
}
 
sno7 <- c()
for (i in row.names(noncoding)) 
  {print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if(a1==b1 & 
     a2==b2)
  { sno7 <-rbind(sno7,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(b1-a2)/(b1-b2+a1-a2), overlap=b1-a2)
                )}
}
#sno1=class A overlapp,sno2=class B,sno3&sno6=classD, sno4&sno5=class C,sno7=class E)
snoa <- sno1
snob <- sno2
sno4
snod<- rbind(sno4,sno5)
snoc <- rbind(sno3,sno6)
snoe <- sno7
snoe <- rbind(snoa[1,],snoe)
snoe <- rbind(snob[1:5,],snoe)
snoe <- rbind(snoc[1:7,],snoe)
snoa <- snoa[-1,]
snob <- snob[-c(1:5),]
snoc <- snoc[-c(1:7),]

# PCC for class A
noncoding$id <- paste0(noncoding$id1,"_",noncoding$id2)
snoa$id <- paste0(snoa$id1,"_",snoa$id2)
## merge data together and only keep id equial from snob
cornoncoda<- merge(snoa,noncoding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(cornoncoda, file="cor_noncodingA.txt", sep="\t", quote=F, row.names=F, col.names=T)
# PCC for class B 
snob$id <- paste0(snob$id1,"_",snob$id2)
cornoncodb<- merge(snob,noncoding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(cornoncodb, file="cor_noncodingB.txt", sep="\t", quote=F, row.names=F, col.names=T)
# PCC for class C 
snoc$id <- paste0(snoc$id1,"_",snoc$id2)
cornoncodc<- merge(snoc,noncoding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(cornoncodc, file="cor_noncodingC.txt", sep="\t", quote=F, row.names=F, col.names=T)
# PCC for class D
snod$id <- paste0(snod$id1,"_",snod$id2)
cornoncodd<- merge(snod,noncoding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(cornoncodd, file="cor_noncodingD.txt", sep="\t", quote=F, row.names=F, col.names=T)
# PCC for class E
snoe$id <- paste0(snoe$id1,"_",snoe$id2)
cornoncode<- merge(snoe,noncoding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(cornoncode, file="cor_noncodingE.txt", sep="\t", quote=F, row.names=F, col.names=T)
# merge all classes
non1 <- rbind(cornoncoda,cornoncodb)
non2 <- rbind(non1,cornoncodc)
non3 <- rbind(non2,cornoncodd)
noncod<- rbind(non3,cornoncode)
write.table(noncod, file="cor_noncoding06292021.txt", sep="\t", quote=F, row.names=F, col.names=T)
