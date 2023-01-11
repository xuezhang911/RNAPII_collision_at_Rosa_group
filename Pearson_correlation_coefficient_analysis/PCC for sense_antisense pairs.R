# generate PCC for non-annotated and annotated sense/antisense gene pairs in 144 Arabidopsis accessions: DOI https://doi.org/10.1038/nature11968
getwd()
setwd("/Users/xung0001/Desktop")
#load counts 
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
ann <- s[,1:8]
slength <- s$V3-s$V2
id <- s$V4
as <- read.table("counts.antisens.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
as$V4 <- id
ann <- as[,1:8]
aslength <- as$V3-as$V2
# load sample accessions
samples <- as.character(read.table("bams.txt", sep=" ", quote="", header=F, stringsAsFactors=F)[1,])
samples <- gsub("\\.bam","",gsub(".*\\/","",samples))
# prepare matrix
s <- s[,-c(1,2,3,4,5,6)]
s <- as.data.frame(lapply(s, as.numeric))
rownames(s) <- id
colnames(s) <- samples
as <- as[,-c(1,2,3,4,5,6)]
as <- as.data.frame(lapply(as, as.numeric))
rownames(as) <- id
colnames(as) <- samples
as <- as.matrix(as)
s <- as.matrix(s)
# filter out some raw counts
s.median <- median(s[rownames(s)=="AT5G10140",])
as.median <- median(as[rownames(as)=="AT5G10140",])
as[rownames(as)=="AT5G10140",]
cor.test(s[rownames(as)=="AT5G10140",], as[rownames(as)=="AT5G10140",])
cor.test(s[rownames(as)=="AT1G69570",], as[rownames(as)=="AT1G69570",])
# colum 69=Je-0
s <- s[,-69]
as[,69:70]
as <- as[,-69]
# correlation anlysis
test_f <- data.frame(NULL)
for(i in 1:nrow(s))
{
  print(i)
  s.mean <- mean(s[i,])
  s.median <- median(s[i,])
  as.mean <- mean(as[i,])
  as.median <- median(as[i,])
  
  if(s.median >=3 & as.median >= 3) {
    ct <- cor.test(s[i,], as[i,])
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
test_f <- test_f[order(test_f$cor, decreasing=T),]
test_s <- test_f[order(test_f$V4, decreasing=F),]
write.table(test_f, file="test_f.txt", sep="\t", quote=F, row.names=F, col.names=T)
test_f <- read.table("test_f.txt", sep="\t", header = T)
Cor <- test_f
## load sense-antisense pairs
co_antlnc <- read.table("coding-antilncoding.txt", sep=",",header = T)
co_coding <- read.table("coding_anticoding.txt", sep=",", header=T)
co_noncoding <- read.table("coding-noncoding.txt", sep=",", header = T)
non_noncoding <- read.table("noncoding_noncoding.txt", sep=",", header = T)

co_antlnc$V2[order(co_antlnc$V2)]
co_coding$V2[order(co_coding$V2)]
# extract coding_antlncRNAantisense information from correlation test
cor_co_antlnc <- data.frame(NULL)
for (i in intersect(co_antlnc$V2, Cor$V4))
{print(i)
  cor_co_antlnc <- rbind(cor_co_antlnc, Cor[Cor$V4==i, ])
}
cor_co_antlnc <- na.omit(cor_co_antlnc[order(cor_co_antlnc$cor,decreasing = T),])
write.table(cor_co_antlnc, file="cor_co-antlncRNA.txt", sep="\t", quote=F, row.names=F, col.names=T)
cor_co_antlncRNA <- read.table("cor_co-antlncRNA.txt", sep="\t",header =T)
# extract coding_noncoding information from correlation test
cor_co_nonco <- data.frame(NULL)
for (i in intersect(co_noncoding$V2, Cor$V4))
{print(i)
  cor_co_nonco <- rbind(cor_co_nonco, Cor[Cor$V4==i, ])
}
cor_co_nonco <- na.omit( cor_co_nonco [order( cor_co_nonco $cor,decreasing = T),])
write.table( cor_co_nonco , file="cor_co_noncoding.txt", sep="\t", quote=F, row.names=F, col.names=T)
cor_co_noncoding <- read.table("cor_co_noncoding.txt", sep="\t",header = T)
# extract coding_coding information from correlation test
cor_co_coding <- data.frame(NULL)
for (i in intersect(co_coding$V2, Cor$V4))
{print(i)
  cor_co_coding <- rbind(cor_co_coding , Cor[Cor$V4==i, ])
}
cor_co_coding <- na.omit(cor_co_coding [order( cor_co_coding$cor,decreasing = T),])
write.table( cor_co_coding  , file="cor_co_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
cor_co_coding <- read.table("cor_co_coding.txt", sep="\t",header = T)
# extract noncoding-noncoding 
cor_non_noncoding <- data.frame(NULL)
for (i in intersect(non_noncoding$V2, Cor$V4))
{print(i)
  cor_non_noncoding  <- rbind(cor_non_noncoding , Cor[Cor$V4==i, ])
}
cor_non_noncoding <- na.omit(cor_non_noncoding [order( cor_non_noncoding$cor,decreasing = T),])
write.table( cor_co_nonco , file=" cor_co_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
cor_co_coding <- cor_co_coding[order(cor_co_coding$p.value,decreasing = F),]

cor_co_coding <- cor_co_coding[which(cor_co_coding$p.value<=0.01),]
cor_co_coding <- cor_co_coding[order(cor_co_coding$p.value,decreasing = F),]
hist(cor_co_coding$p.value)
dev.new()

cor_co_nonco <-cor_co_nonco [order(cor_co_nonco$p.value,decreasing = F),]
cor_co_nonco$p.value<=0.01
cor_co_nonco <-cor_co_nonco[which(cor_co_nonco$p.value<=0.01),]
