getwd()
setwd("/Users/xung0001/Desktop")
## load pair-wise pearson correlation analysis result
ctrco<- read.table("cor_coding_coding.txt", sep="\t",header = T)
ctrco <- ctrco[ctrco$p.value<=0.01,]
ctrco <- ctrco[order(ctrco$cor,decreasing = T),]
ctrco0.6 <- ctrco[which(ctrco$cor>=0.6),]
ctr<- read.table("cor_coding_noncoding.txt", sep="\t",header = T)
ctr <- ctr[ctr$p.value<=0.01,]
ctr <- ctr[order(ctr$cor,decreasing = T),]
#load libraries
library(clusterProfiler)
library(AnnotationHub)
library(ggplot2)
library(stringr)
library(DOSE)
library(pathview)
library(topGO)
library(dplyr)
library(BiocFileCache)
library(org.At.tair.db)
##KEGG 
genelist <-rbind(corcneu$id1,corcneu$id2)
kk <- enrichKEGG(gene=genelist,organism="ath", pAdjustMethod = "BH", pvalueCutoff=0.99,qvalueCutoff=0.99) # you can set cut-off by yourself
kk <- as.data.frame(kk)
colnames(kk)
kk <- kk[order(kk$p.adjust,decreasing = F),]
g_kegg <- ggplot(kk, aes(x=reorder(Description,order(p.adjust, decreasing = T)), y=Count)) + 
  geom_bar(aes(fill=p.adjust),stat="identity", width=0.8,colour="white",position=position_dodge()) +
  scale_x_discrete(name ="Pathway names",breaks=) +
  scale_y_continuous(name ="gene number") +
  coord_flip() +
  ggtitle("Pathway Enrichment")+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),axis.text=element_text(size=12,colour = "black",vjust = 0.5, hjust = 0.5),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"))
print(g_kegg)
head(kk)
 ## save as pdf with the changed settings of width and height
write.csv(kk,file = "kegg_pearsoncoding_noncoding.csv")
##Go
genelist <- ctr$id1
genelist <- nonB$id1
for (onto in c('CC','BP','MF')){
  
  ego <- enrichGO(gene         = genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           =  onto,
                  keyType = "TAIR",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99)}
head(ego)
## separately
egobp <- enrichGO(gene         = genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           =  "BP",
                  keyType = "TAIR",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99)

egocc <- enrichGO(gene         = genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           =  "CC",
                  keyType = "TAIR",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99)
 
egomf <- enrichGO(gene         = genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           =  "MF",
                  keyType = "TAIR",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99)

dotplot(egobp)
dotplot(egocc)
dotplot(egomf)

## ego <- setReadable(ego, org.At.tair.db ) if there is no symbolname of genes, gene name will return NA
write.csv(as.data.frame(ego),file = "goEN.csv")
## go visualization
ego2 <- gofilter(ego,4)
 ego2 <- as.data.frame(ego2)
colnames(ego2)
str_sub(ego2$Description,start = 1,end = 20)
ego2$Description <- str_sub(ego2$Description,start = 1,end = 20)
p <- ggplot(ego2, aes(x=reorder(Description,order(p.adjust, decreasing = T)), y=Count)) + 
  geom_bar(aes(fill=p.adjust),stat="identity", width=0.8,colour="grey",position=position_dodge()) + 
  xlab("Go term") +
  ylab("gene number") +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),axis.text=element_text(size=12,colour = "black",vjust = 0.5, hjust = 0.5),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"))
print(p)
plotGOgraph(egomf)

