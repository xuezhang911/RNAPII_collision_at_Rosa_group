# check the PCC of my candidates of coding-noncoding set in 728 Arabidopsis accessions
setwd("/Users/xung0001/Desktop")
## load normalized counts from 728 accessions except Col0
h <- read.csv("GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv",quote = " ",header=T,sep="\t")
row.names(h) <- h$gene_id
h <- h[,-1]
h <- as.matrix(h)
## load sense-antisense pairs
co_noncoding <- read.table("coding-noncoding.txt", sep=",", header = T)
# select s.median of coding-coding genes greater than 3
s.median <- c()
for (i in row.names(co_noncoding))
{print(i)
  s.median <-rbind(s.median,data.frame(id1=co_noncoding[i,1],id2=co_noncoding[i,2], s.median1=median(h[which(row.names(h)==co_noncoding[i,1]),]),
                                       s.median2=median(h[which(row.names(h)==co_noncoding[i,2]),])))  
}

s.median <- na.omit(s.median)
s.median <-  s.median[which(s.median$s.median1>=3),]
s.median <-  s.median[which(s.median$s.median2>=3),]
### correlation analysis for different groups coding_noncoding
ctr <- data.frame(NULL)
for (i in row.names(s.median))
{print(i)
  ct <- cor.test(h[s.median[i,1],],h[s.median[i,2],])
  mean1 <-median(h[which(row.names(h)==s.median[i,1]),])
  mean2 <- median(h[which(row.names(h)==s.median[i,2]),])
  ctr <- rbind(
    ctr,data.frame(id1=s.median[i,1], id2=s.median[i,2], cor=ct$estimate,p.value=ct$p.value,
                   stringsAsFactors=F, check.names=F ,mean1=mean1,mean2=mean2, ratio=mean2/mean1
    ))}
hist(ctr$p.value,fre=F)
# notice p.value has been adjusted automatically by cor.test
ctr <- ctr[order(ctr$p.value,decreasing = F),]
ctr <- ctr[ctr$p.value<=0.05,]
ctr <- ctr[order(ctr$cor,decreasing = T),]
# another table for correlation test
write.table(ctr, file="cor_noncoding_dec09.txt", sep="\t", quote=F, row.names=F, col.names=T)
