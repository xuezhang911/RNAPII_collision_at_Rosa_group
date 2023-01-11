#### generating coding/antisense coding pairs; coding/antisense non-coding pairs
#library we need
# to install GenomicRanges, it requires S4vectors
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("S4Vectors")
# install GenomicRanges
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicRanges")

library(GenomicRanges)
genes_araport_adj <- readRDS("genes_araport_adj.RDS")
mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]
out1<- na.omit(read.table("out1.txt",header = T,sep = ","))
out2 <- na.omit(read.table("out2.txt",header = T,sep = ","))
out3<- na.omit(read.table("out3.txt",header = T,sep = ","))

out1 <- mRNA[out1$gene_id,]
out2 <- mRNA[out2$gene_id,]
out3 <- mRNA[out3$gene_id,]
out1 <- out1[strand(out1)=="+",]
head(out1)
# Switch the strand of transcripts of out2 or out3
# shift ourt2 genes to 1000bp upstream : for examples the initial arrange is 6788-9130,now 5788-9130
out2 <- shift(out2,-1000) 
out2_sw <- out2
strand(out2_sw) <- ifelse(strand(out2_sw) == "+", "-", "+")
out2_sw_start <- resize(out2_sw, 1, "end")
# Find overlaps with genes in the antisense orientation: for out1 and out2
over_as <- out2_sw %over% out1
head(over_as)
hits <- findOverlaps(out2_sw, out1,type="any")
head(hits)
out2_sw_par <- out2_sw_start[queryHits(hits)]
head(out2_sw_par)
out1_par <- out1[subjectHits(hits)]
coding_antico_1 <- cbind(out1_par$gene_id,out2_sw_par$gene_id)
# after we found all convergent pairs between out1 and out2, we wanted to know how much of them overlaps with antisense lncRNAs
noncoding <- read.table("coding-noncoding.txt", sep=",", quote=",", header=F, stringsAsFactors=F)
noncoding <- noncoding[-6,]
head(noncoding)
coding_antico_1 <- as.data.frame(coding_antico_1)
rt1_1 <- intersect(noncoding$V1,coding_antico_1$V1)
rt1_2 <- intersect(noncoding$V1,coding_antico_1$V2)
#rt1_1+rt1_2=160+157=317
length(intersect(noncoding$V1,coding_antico_1$V1))
#when extract the "-" strand of out1
out1<- na.omit(read.table("out1.txt",header = T,sep = ","))
out1 <- mRNA[out1$gene_id,]
out1 <- out1[strand(out1)=="-",]
head(out1)
out2 <- na.omit(read.table("out2.txt",header = T,sep = ","))
out2 <- mRNA[out2$gene_id,]
head(out2)
# Switch the strand of transcripts of out2 or out3
# shift ourt2 genes to 250 downstream : try to remove ones with only less than 250bp overlap 
out2 <- shift(out2,250)
out2_sw <- out2
head(out2)
strand(out2_sw) <- ifelse(strand(out2_sw) == "+", "-", "+")

head(out2_sw_start)
head(out2_sw)
head(out1)
# Find overlaps with genes in the antisense orientation: for out1 and out2
over_as <- out2_sw %over% out1
head(over_as)
hits <- findOverlaps(out2_sw , out1,type="any")
head(hits)
out2_sw_par <- out2_sw_start[queryHits(hits)]
head(out2_sw_par)
out1_par <- out1[subjectHits(hits)]
out2_sw_start <- resize(out2_sw, 1, "end")
coding_antico_2 <- cbind(out1_par$gene_id,out2_sw_par$gene_id)
coding_antico_2 <- as.data.frame(coding_antico_2)
length(intersect(noncoding$V1,coding_antico_2$V1))
rt2_1<- intersect(noncoding$V1,coding_antico_2$V1)
rt2_1
rt2_2<- intersect(noncoding$V1,coding_antico_2$V2)
rt2_2
#"AT3G27810" "AT1G05730no" "AT2G11851no" "AT2G32650" "AT2G46570" "AT3G24514" "AT5G41860no" "AT5G44920"
#8+8-3=13
## for out2 and out3
# when out2 has + strand
out2<- na.omit(read.table("out2.txt",header = T,sep = ","))
out2 <- mRNA[out2$gene_id,]
out2 <- out2[strand(out2)=="+",]
head(out2)
out3 <- na.omit(read.table("out3.txt",header = T,sep = ","))
out3 <- mRNA[out3$gene_id,]
head(out3)
# shift ourt3 genes to 1000bp upstream : for examples the initial arrange is 6788-9130,now 5788-9130
out3<- shift(out3,-1000) 
head(out3)
out3_sw <- out3
strand(out3_sw) <- ifelse(strand(out3_sw) == "+", "-", "+")
out3_sw_start <- resize(out3_sw, 1, "end")
over_as <- out3_sw %over% out2
hits <- findOverlaps(out3_sw , out2)
out3_sw_par <- out3_sw_start[queryHits(hits)]
out2_par <- out2[subjectHits(hits)]
coding_antico_3 <- cbind(out2_par$gene_id,out3_sw_par$gene_id)
coding_antico_3 <- as.data.frame(coding_antico_3)
length(intersect(noncoding$V1,coding_antico_3$V1))
rt3_1 <- intersect(noncoding$V1,coding_antico_3$V1)
rt3_2 <- intersect(noncoding$V1,coding_antico_3$V2)
head(rt3_2)
#rt3=148+139=287
# when out2 has - strand
out2<- na.omit(read.table("out2.txt",header = T,sep = ","))
out2 <- mRNA[out2$gene_id,]
out2 <- out2[strand(out2)=="-",]
head(out2)
out3 <- na.omit(read.table("out3.txt",header = T,sep = ","))
out3 <- mRNA[out3$gene_id,]
head(out3)
out3 <- shift(out3,250)
out3_sw <- out3
strand(out3_sw) <- ifelse(strand(out3_sw) == "+", "-", "+")
out3_sw_start <- resize(out3_sw, 1, "end")
over_as <- out3_sw %over% out2
hits <- findOverlaps(out3_sw , out2)
out3_sw_par <- out3_sw_start[queryHits(hits)]
out2_par <- out2[subjectHits(hits)]
coding_antico_4 <- cbind(out2_par$gene_id,out3_sw_par$gene_id)
coding_antico_4 <- as.data.frame(coding_antico_4)
length(intersect(noncoding$V1,coding_antico_4$V1))
rt4 <- intersect(noncoding$V1,coding_antico_4$V1)
rt4_2<- intersect(noncoding$V1,coding_antico_4$V2)
rt4_2

#rt4=10_5-2-1=12
#AT2G31900no" "AT3G13530" "AT1G03410no" "AT1G06515" "AT1G68940" "AT2G07692" "AT2G22090" "AT4G00760" "AT4G33300" "AT5G15840"
#"AT2G31910no" "AT1G06500" "AT1G64630" "AT2G16365" "AT3G19590"
totalrt <-length(rt1)+length(rt2)+length(rt3)+length(rt4)
# 635pairs have the arrangment of coding-antisense noncoding-coding out of 1486 pairs
noncoding1 <- noncoding[noncoding$V1%in%rt1_1,]
noncoding1 <- noncoding1[order(noncoding1$V1),]
noncoding1 <- unique(noncoding1)
noncoding2 <- noncoding[noncoding$V1%in%rt1_2,]
noncoding2 <- noncoding2[order(noncoding2$V1),]
noncoding2 <- unique(noncoding2)
noncoding3 <- noncoding[noncoding$V1%in%rt3_1,]
noncoding3 <- noncoding3[order(noncoding3$V1),]
noncoding3 <- unique(noncoding3)
noncoding4 <- noncoding[noncoding$V1%in%rt3_2,]
noncoding4 <- noncoding4[order(noncoding4$V1),]
noncoding4 <- unique(noncoding4)
noncodin1_2 <- noncoding[noncoding$V1%in%rt2_1,]
noncodin1_2 <- noncodin1_2[order(noncodin1_2$V1),]
noncodin1_2 <- unique(noncodin1_2)
noncoding2_2<- noncoding[noncoding$V1%in%rt2_2,]
noncoding2_2 <- noncoding2_2[order(noncoding2_2$V1),]
noncoding2_2 <- unique(noncoding2_2)
noncodin4_1 <- noncoding[noncoding$V1%in%rt4,]
noncodin4_1 <- noncodin4_1[order(noncodin4_1$V1),]
noncodin4_1 <- unique(noncodin4_1)
noncodin4_2 <- noncoding[noncoding$V1%in%rt4_2,]
noncodin4_2 <- noncodin4_2 [order(noncodin4_2$V1),]
noncodin4_2 <- unique(noncodin4_2)
noncoding1_2 <- rbind(noncoding1,noncoding2)
noncoding1_2 <- rbind(noncoding1_2,noncodin1_2)
noncoding1_2 <- rbind(noncoding1_2,noncoding2_2)

noncoding3_4 <- rbind(noncoding3,noncoding4)
noncoding3_4 <- rbind(noncoding3_4,noncodin4_1)
noncoding3_4 <- rbind(noncoding3_4,noncodin4_2)
sense_read_through_antisense <- rbind(noncoding1_2,noncoding3_4)
sense_read_through_antisense <- unique(sense_read_through_antisense[order(sense_read_through_antisense$V1),])
#642pairs have read-through arrangement
write.table(sense_read_through_antisense,"csense_read_through_antisensepair.txt",quote =F,row.names = F,sep = ",")
h <- na.omit(read.table("csense_read_through_antisensepair.txt",header = T,sep = ","))
#distance smaller than 1000bp convergent gene pairs
coding_antico_1_2 <- unique(rbind(coding_antico_1, coding_antico_2))
coding_antico_3_4 <-unique(rbind(coding_antico_3, coding_antico_4))
convergent_coding_pairs <- rbind(coding_antico_1_2,coding_antico_3_4)
# in total there are 3959 pairs of convergent genes (less than 1000bp)
write.table(convergent_coding_pairs,"convergent_coding_less1000_pairs.txt",quote =F,row.names = F,sep = ",")
con <- convergent_coding_pairs
## for coding-noncoding pairs, how many of them have coding-anticoding arrangements
dim(noncoding)
noncoding$convergentdowngene <- c(NA)
for (i in 1:dim(noncoding)[1]) {
  print(i)
  if (noncoding$V1[i]%in%con$V1) {
    noncoding$convergentdowngene[i] <- con[con$V1==noncoding$V1[i],]$V2} 
  if (noncoding$V1[i]%in%con$V2) {noncoding$convergentdowngene[i] <- con[con$V2==noncoding$V1[i],]$V1}}
noncoding <- noncoding[order(noncoding$convergentdowngene),]
nonco <- na.omit(noncoding)
write.table(noncoding,"coding_noncoding_1486_arrangement.txt",quote =F,row.names = F,sep = ",")
## try to select coding-coding overlapping genes from this list
cod <-  read.table("coding_anticoding.txt",header=T,sep = ",")
s <- read.table("counts.sense.txt", sep="\t", quote="", header=T, stringsAsFactors=F)
slength <- s$end-s$start
s <- s[!s$chrom=="chloroplast",]
s <- s[!s$chrom=="mitochondria",]
s<- s[,1:8]
#add gene length
nonco$lengthup2 <- c(NA)
nonco$lengthup1 <- c(NA)
for (i in 1:dim(nonco)[1])
{print(i)
 nonco$lengthup2[i]<-s[which(s$name==nonco[i,3]),]$end-s[which(s$name==nonco[i,3]),]$start
 nonco$lengthup1[i]<-s[which(s$name==nonco[i,1]),]$end-s[which(s$name==nonco[i,1]),]$start
}

#add overlap length between convergent genes
nonco$overlap <- c(NA)
for (i in 1:dim(nonco)[1])
{ print(i)
  a1<- s[s$name==nonco[i,3],]$start 
  a2 <-  s[s$name==nonco[i,3],]$end 
  b1 <-  s[s$name==nonco[i,1],]$start
  b2 <- s[s$name==nonco[i,1],]$end 
  if (a2<=b2 & a1>=b1 & a2>b1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="A",overlap=a2-a1))}
  if (a2>=b2 & a1<=b1 & b1<a2) {overlap<- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="B",overlap=b2-b1))}
  if (a2>b2 & a1<b2 & b1<a1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="C",overlap=b2-a1))}
  if (a2<b2 & a1<b1 & b1<a2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="D",overlap=a2-b1))}
}
## noncoding/non-coding gene pairs
# load the file of coding-noncoding
nonc <- read.table("cor_noncoding_updata.txt",header=T)
dim(nonc)
nonc$convergentdowngene <- c(NA)
for (i in 1:dim(nonc)[1]) {
  print(i)
  if (nonc$id1[i]%in%con$V1) {
    nonc$convergentdowngene[i] <- con[con$V1==nonc$id1[i],]$V2} 
  if (nonc$id1[i]%in%con$V2) {nonc$convergentdowngene[i] <- con[con$V2==nonc$id1[i],]$V1}}
# how many genes of my candidates have the read-through arrangements (24 candidates)
dim(nonc[-which(is.na(nonc$convergentdowngene)),])
y <- nonc[-which(is.na(nonc$convergentdowngene)),][,15]
y
#"AT2G01400" "AT2G45540" "AT1G64570" "AT1G68862" "AT2G20495" "AT4G28640" "AT1G78260" "AT3G02920" "AT2G45360" "AT2G46567"
#"AT3G15200" "AT3G48131" "AT3G63450" "AT3G62080" "AT1G79710" "AT5G19210" "AT4G38260" "AT5G01170" "AT2G31910" "AT4G26500"
# "AT1G72060" "AT1G17870" "AT1G28690" "AT5G59720"
nonc <- nonc[order(nonc$convergentdowngene),]
write.table(nonc,"coding_noncoding_update.txt",quote =F,row.names = F,sep = ",")

nonj <- read.table("coding_noncoding_update.txt",sep=",",header=T)
