##generate correlation value for sense/antisense groups including (coding/coding, coding/noncoding,noncoding/noncoding, linterRNA/mRNA)  and sense/sense overlapping coding/noncoding pairs
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
co_coding <- read.table("coding_anticoding.txt", sep=",", header=T)
# select s.median of coding-coding genes greater than 3
s.median <- c()
for (i in row.names(co_coding))
{print(i)
  s.median <-rbind(s.median,data.frame(id1=co_coding[i,1],id2=co_coding[i,2], s.median1=median(s[which(row.names(s)==co_coding[i,1]),]),
                                       s.median2=median(s[which(row.names(s)==co_coding[i,2]),])))  
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
# set cut_off to keep data 0.05/number of test=784
ctr <- ctr[ctr$p.value<=0.05/dim(ctr)[1],]
ctr <- ctr[order(ctr$cor,decreasing = T),] 
#591 pairs
## extract overalpping information for coding-coding 
s1 <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
s1 <- s1[!s1$chrom=="chloroplast",]
s1 <- s1[,-69]
row.names(s1) <- s1$name
coding <- ctr
sno <- c()
for (i in row.names( coding))
{print(i)
  a1 <-  s1[row.names(s1)==coding[i,2],]$end 
  b1 <-  s1[row.names(s1)==coding[i,1],]$end 
  a2 <- s1[row.names(s1)==coding[i,2],]$start 
  b2 <- s1[row.names(s1)==coding[i,1],]$start
  if( a1>b1 & 
      a2>b2 )
  { sno <-rbind(sno,data.frame(id1=coding[i,1],id2=coding[i,2], overlapp=2*(b1-a2)/(b1-b2+a1-a2), overp=b1-a2)
  )}
}

 sno1<- c()
for (i in row.names(coding))
{print(i)
  a1 <-  s1[row.names(s1)==coding[i,2],]$end 
  b1 <-  s1[row.names(s1)==coding[i,1],]$end 
  a2 <- s1[row.names(s1)==coding[i,2],]$start 
  b2 <- s1[row.names(s1)==coding[i,1],]$start
  if( b1>a1 & 
      b2>a2 )
  { sno1 <-rbind(sno1,data.frame(id1=coding[i,1],id2=coding[i,2], overlapp=2*(a1-b2)/(b1-b2+a1-a2), overp=a1-b2)
  )}
}

sno2<- c()
for (i in row.names(coding))
{print(i)
  a1 <-  s1[row.names(s1)==coding[i,2],]$end 
  b1 <-  s1[row.names(s1)==coding[i,1],]$end 
  a2 <- s1[row.names(s1)==coding[i,2],]$start 
  b2 <- s1[row.names(s1)==coding[i,1],]$start
  if( b1>a1 & 
      b2<a2 )
  { sno2 <-rbind(sno2,data.frame(id1=coding[i,1],id2=coding[i,2], overlapp=(a1-a2)/(b1-b2), overp=a1-a2)
  )}
}
sno3<- c()
for (i in row.names(coding))
{print(i)
  a1 <-  s1[row.names(s1)==coding[i,2],]$end 
  b1 <-  s1[row.names(s1)==coding[i,1],]$end 
  a2 <- s1[row.names(s1)==coding[i,2],]$start 
  b2 <- s1[row.names(s1)==coding[i,1],]$start
  if( a1>=b1 & 
      b2>=a2 )
  { sno3 <-rbind(sno3,data.frame(id1=coding[i,1],id2=coding[i,2], overlapp=(b1-b2)/(a1-a2), overp=b1-b2)
  )}
}

snob <- unique(rbind(sno,sno1))

## pcc for class B 
#creat a new variable for corcoding and snob
coding$id <- paste0(coding$id1,"_",coding$id2)
snob$id <- paste0(snob$id1,"_",snob$id2)
## merge data together and only keep id equial from snob
corcodB<- merge(snob,coding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(corcodB, file="cor_codingB.txt", sep="\t", quote=F, row.names=F, col.names=T)
# class A 
snoa <- unique(rbind(sno2,sno3))
snoa$id <- paste0(snoa$id1,"_",snoa$id2)
corcodA <- merge(snoa,coding[,-c(1:2)],by="id",all.x = T)[,-1]
write.table(corcod, file="cor_codingA.txt", sep="\t", quote=F, row.names=F, col.names=T)
ctr <- rbind(corcodA,corcodB)
ctr <- ctr[order(ctr$cor,decreasing = T),] 
write.table(ctr, file="cor_codingcoding06292021.txt", sep="\t", quote=F, row.names=F, col.names=T)
