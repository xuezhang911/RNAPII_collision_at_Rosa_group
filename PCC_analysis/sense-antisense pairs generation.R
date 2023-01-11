
#### generating coding/antisense coding pairs; coding/antisense non-coding pairs
setwd("/Users/xung0001/Desktop")
options(stringsAsFactors=F)
#library we need
library(GenomicRanges)
# open adjusted araport11 Large Granges from desktop from sebastian group 
genes_araport_adj <- readRDS("genes_araport_adj.RDS")
unique(factor(genes_araport_adj$tx_type))
# extract subgroup Granges and their IDs from  genes_araport_adj
senslncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="lnc_RNA"),]
antilncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="antisense_lncRNA"),]
lnc_rna <- c(senslncRNA,antilncRNA)
mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]
tRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="tRNA"),]
ncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="ncRNA"),]
rRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="rRNA"),]
snoRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="snoRNA"),]
snRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="snRNA"),]
noncodingRNA <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA,lnc_rna)
noncodingRNA1 <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA,antilncRNA)
noncodingRNA2 <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA)
noncodingRNA3 <- c(ncRNA,rRNA,snoRNA,snRNA)
noncodingRNA4 <- c(rRNA,snoRNA,snRNA)
noncodingRNA5 <- c(rRNA,snRNA)

##Find  overlaps between mRNA and antisense lncRNA  (either transcript or gene is covered by at least
#Switch the strand of transcripts:
lnc_sw <- lnc_rna
strand(lnc_sw) <- ifelse(strand(lnc_rna) == "+", "-", "+")
lnc_sw_start <- resize(lnc_sw, 1, "end")
# Find overlaps with genes in the antisense orientation:
over_as <- lnc_sw %over% mRNA
hits2 <- disjoin(lnc_sw,mRNA)

hits2 <- findOverlaps(lnc_sw, mRNA)
lnc_sw_start_par <- lnc_sw_start[queryHits(hits2)]
mRNA_par2 <- mRNA[subjectHits(hits2)]
# coding/antisense non-coding  pairs
cod_anlnc <- cbind(mRNA_par2$gene_id, lnc_sw_start_par$gene_id)
write.table(cod_anlnc,"coding-lncRNA.txt",quote =F,row.names = F,sep = ",")
## ##Find  overlaps between mRNA and noncoding at the antisense strand
non_sw <- noncodingRNA
strand(non_sw) <- ifelse(strand(non_sw) == "+", "-", "+")
non_sw_start <- resize(non_sw, 1, "end")
over_as <- non_sw %over% mRNA
hits2 <- findOverlaps(non_sw, mRNA)
non_sw_start_par <- non_sw_start[queryHits(hits2)]
mRNA_par2 <- mRNA[subjectHits(hits2)]
#coding/non-coding gene pairs
cod_noncoding <- cbind(mRNA_par2$gene_id, non_sw_start_par$gene_id)
write.table(cod_noncoding,"coding-noncoding.txt",quote =F,row.names = F,sep = ",")
## find overlaps btween sense coding and antisense coding genes
# generating three mRNAs subgroups
mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[order(mRNA$gene_id,decreasing=F),]
out1 <- c()
out2 <- c()
out3 <- c()
for (i in c(0:31556))
{ print(i)
  out1 <- rbind(out1,mRNA[1+3*i,])
  out2 <-rbind(out2,mRNA[2+3*i,])
  out3 <- rbind(out3,mRNA[3+3*i,])
}
write.table(out1,"out1.txt",quote =F,row.names = F,sep = ",")
write.table(out2,"out2.txt",quote =F,row.names = F,sep = ",")
write.table(out3,"out3.txt",quote =F,row.names = F,sep = ",")
out1<- na.omit(read.table("out1.txt",header = T,sep = ","))
out2 <- na.omit(read.table("out2.txt",header = T,sep = ","))
out3<- na.omit(read.table("out3.txt",header = T,sep = ","))
## generating coding/antisense coding genes
# extract out1,out2,out3 granges objects from mRNA
mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]
out1 <- mRNA[out1$gene_id,]
out2 <- mRNA[out2$gene_id,]
out3 <- mRNA[out3$gene_id,]
# Switch the strand of transcripts of out2 or out3
out2_sw <- out2
strand(out2_sw) <- ifelse(strand(out2_sw) == "+", "-", "+")
out2_sw_start <- resize(out2_sw, 1, "end")
# Find overlaps with genes in the antisense orientation: for out1 and out2
over_as <- out2_sw %over% out1
hits <- findOverlaps(out2_sw , out1)
out2_sw_par <- out2_sw_start[queryHits(hits)]
out1_par <- out1[subjectHits(hits)]
coding_antico_1 <- cbind(out1_par$gene_id,out2_sw_par$gene_id)
# for out2 and out3
out3_sw <- out3
strand(out3_sw) <- ifelse(strand(out3_sw) == "+", "-", "+")
out3_sw_start <- resize(out3_sw, 1, "end")
over_as <- out3_sw %over% out2
hits <- findOverlaps(out3_sw , out2)
out3_sw_par <- out3_sw_start[queryHits(hits)]
out2_par <- out2[subjectHits(hits)]
coding_antico_2 <- cbind(out2_par$gene_id,out3_sw_par$gene_id)
# coding_anti coding mRNA pairs
coding_antico <- rbind(coding_antico_1,coding_antico_2)
write.table(coding_antico,"coding_anticoding.txt",quote =F,row.names = F,sep = ",")
## noncoding/non-coding gene pairs
#  Switch the strand of transcripts of noncodingRNA1 overlap with senselncRNA
nc_sw <-noncodingRNA1
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% senslncRNA
hits <- findOverlaps(nc_sw , senslncRNA)
nc_sw_par <- nc_sw_start[queryHits(hits)]
osen_par <- senslncRNA[subjectHits(hits)]
nonc_annonc1 <- cbind(osen_par$gene_id,nc_sw_par$gene_id)
#  Switch the strand of transcripts of noncodingRNA2 overlap with antilncRNA
nc_sw <-noncodingRNA2
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% antilncRNA
hits <- findOverlaps(nc_sw , antilncRNA)
anti_sw_par <- nc_sw_start[queryHits(hits)]
sen_par <- antilncRNA[subjectHits(hits)]
nonc_annonc2 <- cbind(sen_par$gene_id,anti_sw_par$gene_id)
# overlap noncoding3 with tRNA #0 result
nc_sw <-noncodingRNA3
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% tRNA
hits <- findOverlaps(nc_sw , tRNA)
anti_sw_par <- nc_sw_start[queryHits(hits)]
tRNA_par <- tRNA[subjectHits(hits)]
nonc_annonc3 <- cbind(tRNA_par$gene_id,anti_sw_par$gene_id)
# overlap noncoding4 with ncRNA #0 result
nc_sw <-noncodingRNA4
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% ncRNA
hits <- findOverlaps(nc_sw , ncRNA)
anti_sw_par <- nc_sw_start[queryHits(hits)]
ncRNA_par <- ncRNA[subjectHits(hits)]
nonc_annonc34<- cbind(ncRNA_par$gene_id,anti_sw_par$gene_id)
# overlap noncoding5 with snoRNA #0 result
nc_sw <-noncodingRNA5
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% snoRNA
hits <- findOverlaps(nc_sw , snoRNA)
anti_sw_par <- nc_sw_start[queryHits(hits)]
snoRNA_par <- snoRNA[subjectHits(hits)]
nonc_annonc5<- cbind(snoRNA_par$gene_id,anti_sw_par$gene_id)
# overlap rRNA with snRNA #0 result
nc_sw <-snRNA
strand(nc_sw)<- ifelse(strand(nc_sw) == "+", "-", "+")
nc_sw_start <- resize(nc_sw, 1, "end")
over_as <- nc_sw %over% rRNA
hits <- findOverlaps(nc_sw , rRNA)
anti_sw_par <- nc_sw_start[queryHits(hits)]
rRNA_par <- rRNA[subjectHits(hits)]
nonc_annonc6<- cbind(rRNA_par$gene_id,anti_sw_par$gene_id)
# merge noncoding-antisense noncoding pairs
non_antinoncod <- rbind(cod_noncoding,nonc_annonc1,nonc_annonc2)
write.table(non_antinoncod,"noncoding_noncoding.txt",quote =F,row.names = F,sep = ",")
## sense_sense overlapping gene pairs
# Find overlaps with genes in the sense orientation(sense_sense tandem): for out1 and out2 
over_as <- out2 %over% out1
hits <- findOverlaps(out2, out1)
out2_sw_par <- out2[queryHits(hits)]
out1_par <- out1[subjectHits(hits)]
coding_codingtandem_1 <- cbind(out1_par$gene_id,out2_sw_par$gene_id)
# Find overlaps with genes in the sense orientation(sense_sense tandem): for out2 and out3
over_as <- out3 %over% out2
hits <- findOverlaps(out3, out2)
out3_sw_par <- out3[queryHits(hits)]
out2_par <- out2[subjectHits(hits)]
coding_codingtandem_2 <- cbind(out2_par$gene_id,out3_sw_par$gene_id)
coding_codingtandem <- rbind(coding_codingtandem_1,coding_codingtandem_2)
write.table(coding_codingtandem,"coding_codingtandem.txt",quote =F,row.names = F,sep = ",")
# find overlaps with genes in the sense noncoding
over_as <- mRNA %over% noncodingRNA
hits <- findOverlaps(mRNA, noncodingRNA)
mRNA_sw_par <- mRNA[queryHits(hits)]
noncoding_par <- noncodingRNA[subjectHits(hits)]
coding_noncodingtandem <- cbind(mRNA_sw_par$gene_id,noncoding_par$gene_id)
write.table(coding_noncodingtandem,"coding_noncodingtandem.txt",quote =F,row.names = F,sep = ",")
