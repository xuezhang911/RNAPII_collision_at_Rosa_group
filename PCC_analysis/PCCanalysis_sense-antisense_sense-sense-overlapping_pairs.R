##generate correlation value for sense/antisense groups including (coding/coding, coding/noncoding,noncoding/noncoding, linterRNA/mRNA)  and sense/sense overlapping coding/noncoding pairs
getwd()
setwd("/Users/xung0001/Desktop")
# load raw counts data
#load counts 
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)

s <- s[!s$chrom=="chloroplast",]
row.names(s) <- s$name
ann <- s[,1:8]
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
as_scale<- t(scale(t(as)))
as_scale <- as.matrix(as_scale)
#correlation nalysis
test_f <- data.frame(NULL)
for(i in 1:nrow(s))
{
  print(i)
  s.mean <- mean(s[i,])
  s.median <- median(s[i,])
  as.mean <- mean(as[i,])
  as.median <- median(as[i,])
  
  if(s.median >=3 & as.median >= 3) {
    ct <- cor.test(s_scale[i,], as_scale[i,])
    corr <- ct$estimate
    p.value <- ct$p.value
  }  else {
    corr <- NA
    p.value <- NA
  }
  test_f <- rbind(
    test_f,
    data.frame(ann[i,],
               cor=corr,
               p.value=p.value,
               sens.median=s.median,
               antisens.median=as.median,
               t(s[i,]/as[i,]),
               stringsAsFactors=F, check.names=F)
  )
  
}
slength <- s$end-s$start
rpkm<- function(counts,genelength){
  rate <- counts/genelength
  t(t(rate)/colSums(counts))*10^9
}
sr<- rpkm(s,slength)
sr <- na.omit(sr)
##generate correlation value for sense/antisense groups including (coding/coding, coding/noncoding,noncoding/noncoding, linterRNA/mRNA)  and sense/sense overlapping coding/noncoding pairs
# load coding/coding, coding/noncoding,noncoding/noncoding
## load sense-antisense pairs
co_antlnc <- read.table("coding-antilncoding.txt", sep=",",header = T)
co_coding <- read.table("coding_anticoding.txt", sep=",", header=T)
co_noncoding <- read.table("coding-noncoding.txt", sep=",", header = T)
non_noncoding <- read.table("noncoding_noncoding.txt", sep=",", header = T)
## extract s.median value from co_coding pairs
s.median <- c()
for (i in row.names(co_coding))
{print(i)
  s.median <-rbind(s.median,data.frame(id1=co_coding[i,1],id2=co_coding[i,2], s.median1=median(s[which(row.names(s)==co_coding[i,1]),]),
                                       s.median2=median(s[which(row.names(s)==co_coding[i,2]),])))  
}

s.median <- na.omit(s.median)
s.median <-  s.median[which(s.median$s.median1>=3),]
s.median <-  s.median[which(s.median$s.median2>=3),]
write.table(s.median, file="smedian_co_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
##extract s.median value from co_noncoding pairs
s.median_non <- c()
for (i in row.names(co_noncoding))
{print(i)
  s.median_non <-rbind(s.median_non,data.frame(id1=co_noncoding[i,1],id2=co_noncoding[i,2], 
                                               s.median1=median(s[which(row.names(s)==co_noncoding[i,1]),]),
                                               s.median2=median(s[which(row.names(s)==co_noncoding[i,2]),])))
}
s.median_non<- na.omit( s.median_non)
s.median_non <-   s.median_non[which( s.median_non$s.median1>=3),]
s.median_non <-   s.median_non[which( s.median_non$s.median2>=3),]
write.table(s.median_non, file="smedian_co_noncoding.txt", sep="\t", quote=F, row.names=F, col.names=T)
##extract s.median value from non_noncoding pairs
s.mediann_non <- c()
for (i in row.names(non_noncoding))
{print(i)
  s.mediann_non <-rbind(s.mediann_non,data.frame(id1=non_noncoding[i,1],id2=non_noncoding[i,2], s.median1=median(s[which(row.names(s)==non_noncoding[i,1]),]),
                                                 s.median2=median(s[which(row.names(s)==non_noncoding[i,2]),])))
}
s.mediann_non<- na.omit( s.mediann_non)
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median1>=3),]
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median2>=3),]
## correlation analysis for different groups coding_noncoding
ctr <- data.frame(NULL)
for (i in row.names(s.median_non))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==s.median_non[i,1]),],s_scale[which(row.names(s_scale)==s.median_non[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=s.median_non[i,1], id2=s.median_non[i,2], cor=ct$estimate,p.value=ct$p.value,
                   
                   stringsAsFactors=F, check.names=F )
  )}

ctr <- ctr[order(ctr$p.value,decreasing = F),]
ctr <- ctr[ctr$p.value<=0.01,]
ctr <- ctr[order(ctr$cor,decreasing = T),]
write.table(ctr, file="cor_coding_noncoding.txt", sep="\t", quote=F, row.names=F, col.names=T)
## PCC analysis of lincRNA_mRNA
# loading sense-antisense pairs which has no overlap, long intergenic noncoding RNA
linc <- read.table("lincRNA_mRNA.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
s.mediann_non <- c()
for (i in 1:dim(linc)[1])
{print(i)
  s.mediann_non <-rbind(s.mediann_non,data.frame(id1=linc[i,1],id2=linc[i,2], s.median1=median(s[which(row.names(s)==linc[i,1]),]),
                                                 s.median2=median(s[which(row.names(s)==linc[i,2]),])))
}

s.mediann_non<- na.omit( s.mediann_non)
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median1>=3),]
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median2>=3),]
## correlation analysis for different groups coding_noncoding
ctr <- data.frame(NULL)
for (i in row.names(s.mediann_non))
{print(i)
  ct <- cor.test(s[which(row.names(s)==s.mediann_non[i,1]),],s[which(row.names(s)==s.mediann_non[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=s.mediann_non[i,1], id2=s.mediann_non[i,2], cor=ct$estimate,p.value=ct$p.value,
                   
                   stringsAsFactors=F, check.names=F )
  )}

ctr <- ctr[order(ctr$p.value,decreasing = F),]
ctr <- ctr[ctr$p.value<=0.01,]
ctr <- ctr[order(ctr$cor,decreasing = T),]
# variance for ctr group
var1 <- c()
for (i in 1:nrow(ctr)) {
  var1 <- rbind(var1,var(s[ctr[i,1],]))
  
}
var2 <- c()
for (i in 1:nrow(ctr)) {
  var2 <- rbind(var2,var(s[ctr[i,2],]))
  
}
var <- cbind(var1,var2)
ctr <- cbind(ctr,var)
colnames(ctr)[5] <- "variance1"
colnames(ctr)[6] <- "variance2"
write.table(ctr, file="cor_lincRNA_antisensemRNA.txt", sep="\t", quote=F, row.names=F, col.names=T)
## PCC analysis of sense-sense overlapping coding/noncoding pairs
glncoding<- read.table("geniclncRNA_mRNA.txt", sep=",", quote="", header=T, stringsAsFactors=F)
st<- read.table("downlncRNA_mRNA.txt", sep=",", quote="", header=T, stringsAsFactors=F)
sy<- read.table("uplncRNA_mRNA.txt", sep=",", quote="", header=T, stringsAsFactors=F)
names(sy) <- names(st) <- names(glncoding)
udlnc <- rbind(st,sy,glncoding)
s.mediann_non <- c()
for (i in row.names(udlnc))
{print(i)
  s.mediann_non <-rbind(s.mediann_non,data.frame(id1=udlnc[i,1],id2=udlnc[i,2], s.median1=median(s[which(row.names(s)==udlnc[i,1]),]),
                                                 s.median2=median(s[which(row.names(s)==udlnc[i,2]),])))
}
s.mediann_non<- unique(na.omit( s.mediann_non))
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median1>=3),]
s.mediann_non <-   s.mediann_non[which( s.mediann_non$s.median2>=3),]
s.mediann_non <- unique(s.mediann_non)
## correlation analysis for different groups coding_noncoding
ctr <- data.frame(NULL)
for (i in row.names(s.mediann_non))
{print(i)
  ct <- cor.test(s[which(row.names(s)==s.mediann_non[i,1]),],s[which(row.names(s)==s.mediann_non[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=s.mediann_non[i,1], id2=s.mediann_non[i,2], cor=ct$estimate,p.value=ct$p.value,
                   
                   stringsAsFactors=F, check.names=F )
  )}

ctr <- ctr[order(ctr$p.value,decreasing = F),]
ctr <- ctr[ctr$p.value<=0.01,]
ctr <- ctr[order(ctr$cor,decreasing = T),]
# variance for ctr group
var1 <- c()
for (i in 1:nrow(ctr)) {
  var1 <- rbind(var1,var(s[ctr[i,1],]))
  
}
var2 <- c()
for (i in 1:nrow(ctr)) {
  var2 <- rbind(var2,var(s[ctr[i,2],]))
  
}
var <- cbind(var1,var2)
ctr <- cbind(ctr,var)
colnames(ctr)[5] <- "variance1"
colnames(ctr)[6] <- "variance2"
write.table(ctr, file="cor_sense_senseoverlap.txt", sep="\t", quote=F, row.names=F, col.names=T)

## PCC analysis based on the overlapping percentage beween sense and antisense coding/noncoding groups 
coding <- read.table("cor_coding.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
coding <- coding[,-(5:6)]
noncoding <- read.table("cor_coding_noncoding.txt",sep="\t", quote="", header=T, stringsAsFactors=F)
## extract overalpping information for coding-coding (also coding_noncoding)
## extract overalpping information for coding-coding 
sno <- c()
for (i in row.names(coding))
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
  if( b1>=a1 & 
      b2>=a2 )
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
  if( a1>b1 & 
      b2>a2 )
  { sno3 <-rbind(sno3,data.frame(id1=coding[i,1],id2=coding[i,2], overlapp=(b1-b2)/(a1-a2), overp=b1-b2)
  )}
}

##extract overlapping information for coding-noncodingset1
sno <- c()
for (i in row.names(noncoding))
{print(i)
 a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
 b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1>b1 & 
      a2>b2 )
     { sno <-rbind(sno,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(b1-a2)/(b1-b2+a1-a2), overp=b1-a2)
  )}
}
##set2
sno1 <- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1<b1 & 
      a2>b2)
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
  if( (a2<=b2 & a1>=b1 ))
  { sno2 <-rbind(sno2,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=(b1-b2)/(a1-a2),overlap=b1-b2)
  )}
}
sno4 <- c()
for (i in row.names(noncoding))
{print(i)
  a1 <-  s1[row.names(s1)==noncoding[i,2],]$end 
  b1 <-  s1[row.names(s1)==noncoding[i,1],]$end 
  a2 <- s1[row.names(s1)==noncoding[i,2],]$start 
  b2 <- s1[row.names(s1)==noncoding[i,1],]$start
  if( a1<b1 & 
      a2<b2 )
  { sno4 <-rbind(sno4,data.frame(id1=noncoding[i,1],id2=noncoding[i,2], overlapp=2*(a1-b2)/(b1-b2+a1-a2), overlap=a1-b2)
  )}
}
write.table(sco, file="coding-coding_sequence_overlap.txt", sep="\t", quote=F, row.names=F, col.names=T)
coding$overlap <- sco$overlap
write.table(coding, file="cor_coding_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)

s_scale<- t(scale(t(s)))
s_scale <- as.matrix(s_scale)
ctr <- data.frame(NULL)
for (i in row.names(sno))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==sno[i,1]),],s_scale[which(row.names(s_scale)==sno[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=sno[i,1], id2=sno[i,2], overlapp=sno[i,3],overlap=sno[i,4],cor=ct$estimate,p.value=ct$p.value,
                   mean1=mean(sr[which(row.names(sr)==sno[i,1]),]),
                   mean2=mean(sr[which(row.names(sr)==sno[i,2]),]),
                   stringsAsFactors=F, check.names=F )
  )}
cor_noncoding1 <- ctr 
write.table(ctr, file="cor_noncoding1.txt", sep="\t", quote=F, row.names=F, col.names=T)
ctr <- data.frame(NULL)
for (i in row.names(sno1))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==sno1[i,1]),],s_scale[which(row.names(s_scale)==sno1[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=sno1[i,1], id2=sno1[i,2], overlapp=sno1[i,3],overlap=sno1[i,4],cor=ct$estimate,p.value=ct$p.value,
                   mean1=mean(sr[which(row.names(sr)==sno1[i,1]),]),
                   mean2=mean(sr[which(row.names(sr)==sno1[i,2]),]),
                   
                   stringsAsFactors=F, check.names=F )
  )}
cor_noncoding2 <- ctr 
write.table(ctr, file="cor_noncoding2.txt", sep="\t", quote=F, row.names=F, col.names=T)
ctr <- data.frame(NULL)
for (i in row.names(sno2))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==sno2[i,1]),],s_scale[which(row.names(s_scale)==sno2[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=sno2[i,1], id2=sno2[i,2], overlapp=sno2[i,3],overlap=sno2[i,4],
                   cor=ct$estimate,p.value=ct$p.value,
                   mean1=mean(sr[which(row.names(sr)==sno2[i,1]),]),
                   mean2=mean(sr[which(row.names(sr)==sno2[i,2]),]),
                   stringsAsFactors=F, check.names=F )
  )}
cor_noncoding3 <- ctr 
write.table(ctr, file="cor_noncoding3.txt", sep="\t", quote=F, row.names=F, col.names=T)
ctr <- data.frame(NULL)
for (i in row.names(sno4))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==sno4[i,1]),],s_scale[which(row.names(s_scale)==sno4[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=sno4[i,1], id2=sno4[i,2], overlapp=sno4[i,3],overlap=sno4[i,4],cor=ct$estimate,p.value=ct$p.value,mean1=mean(sr[which(row.names(sr)==sno4[i,1]),]),
                   mean2=mean(sr[which(row.names(sr)==sno4[i,2]),]),
                   
                   stringsAsFactors=F, check.names=F )
  )}
cor_noncoding4 <- ctr 
write.table(ctr, file="cor_noncoding4.txt", sep="\t", quote=F, row.names=F, col.names=T)
cor_noncoding <- rbind(cor_noncoding1,cor_noncoding2,cor_noncoding3,cor_noncoding4)
cor_noncoding <- cor_noncoding[order(cor_noncoding$cor,decreasing = T),]
cor_noncoding <- unique(cor_noncoding)
write.table(cor_noncoding, file="cor_noncoding.txt", sep="\t", quote=F, row.names=F, col.names=T)
## for coding-coding
ctr <- data.frame(NULL)
for (i in row.names(sno0))
{print(i)
  ct <- cor.test(s_scale[which(row.names(s_scale)==sno0[i,1]),],s_scale[which(row.names(s_scale)==sno0[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=sno0[i,1], id2=sno0[i,2], overlap=sno0[i,3],cor=ct$estimate,p.value=ct$p.value,mean1=mean(sr[which(row.names(sr)==sno0[i,1]),]),
                   mean2=mean(sr[which(row.names(sr)==sno0[i,2]),]),
                   
                   stringsAsFactors=F, check.names=F )
  )}
cor_coding <- ctr 
cor_coding <- cor_coding[order(cor_coding$cor,decreasing = T),]
cor_coding <- unique(cor_coding)
write.table(cor_coding, file="cor_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
FO1 <- cor_noncoding1[which(cor_noncoding1$overlapp>=0.95),]
FO2 <- cor_noncoding2[which(cor_noncoding2$overlapp>=0.95),]
FO3 <- cor_noncoding3[which(cor_noncoding3$overlapp>=0.95),]
FO4 <- cor_noncoding4[which(cor_noncoding4$overlapp>=0.95),]
FO <- rbind(FO1,FO2,FO3,FO4)
cor_noncoding01 <- cor_noncoding1[-which(cor_noncoding1$overlapp>=0.95),]
cor_noncoding02<- cor_noncoding2[which(cor_noncoding2$overlapp<=0.95),]
cor_noncoding03 <- cor_noncoding3[-which(cor_noncoding3$overlapp>=0.95),]
cor_noncoding04 <- cor_noncoding4[-which(cor_noncoding4$overlapp>=0.95),]
FO <- FO[order(FO$cor,decreasing = T),]
FO <- unique(FO)
co <- rbind(cor_noncoding02,cor_noncoding03,FO)
co <- co[order(co$cor,decreasing = T),]
barplot(co$cor)
dim(co[which(co$cor>=0.5),])[1]/length(co$id1) ##0.41
dim(co[which(co$cor>=0.6 ),])[1]/length(co$id1)
barplot(cor_noncoding01$cor)
dim(cor_noncoding01[which(cor_noncoding01$cor>=0.5 ),])[1]/length(cor_noncoding01$id1)
corcoding <- read.table("cor_codin.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
cornoncoding <- read.table("cor_noncoding.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
cor <- rbind(corcoding,cornoncoding)
cor <- cor[order(cor$cor,decreasing = T),]


