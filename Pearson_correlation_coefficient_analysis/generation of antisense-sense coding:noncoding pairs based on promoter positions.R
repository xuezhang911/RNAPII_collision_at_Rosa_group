## generate antisense-sense pairs based on promoter positions
setwd("/Users/xung0001/Desktop")
options(stringsAsFactors=F)
#library we need
library(rtracklayer)
library(GenomicRanges)
# open adjusted araport11 Large Granges from desktop from sebastian group 
genes_araport_adj <- readRDS("genes_araport_adj.RDS")
unique(factor(genes_araport_adj$tx_type))
# extract lncRNAs/antilncRNA/mRNA Granges and their IDs from  genes_araport_adj
senslncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="lnc_RNA"),]
antilncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="antisense_lncRNA"),]
lnc_rna <- c(senslncRNA,antilncRNA)
mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]

# Find strong overlaps on the sense strand (either transcript or gene is covered by at least 50% width):
over_s <- lnc_rna %over% mRNA
hits <- findOverlaps(lnc_rna, mRNA)
# find overlapping lncRNA and mRNA transcripts ID
lnc_par <- lnc_rna[queryHits(hits)]
lnc_par$gene_id
mRNA_par <- mRNA[subjectHits(hits)]
mRNA_par$gene_id
# extract strong overlapping == genic sense lncRNA and mRNA pairs
overlap <- width(pintersect(lnc_par, mRNA_par))
over_lnc<- overlap / width(lnc_par)
over_mRNA <- overlap / width(mRNA_par)
strong_overlap <- over_lnc >= 0.5 | over_mRNA >= 0.5
strongmRNA <- mRNA_par[strong_overlap]
strongmRNA$gene_id
lncstrong <- lnc_par[strong_overlap]
lncstrong$gene_id
# strong overlapping antisense -sense overlapping pairs
genicID <- cbind(lncstrong$gene_id,strongmRNA$gene_id)
colnames(genicID) <- c("genic_lncRNA","mRNA")
write.table(genicID,"geniclncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")

# Detect weak upstream and downstream overlaps on the same strand:
weak_overlap <- !strong_overlap
a <- start(lnc_par) <= start(mRNA_par) & end(lnc_par) < end(mRNA_par)
b <- start(lnc_par) > start(mRNA_par) & end(lnc_par) >= end(mRNA_par)
over_up <- ifelse(strand(lnc_par) == "+", a, b)
over_down <- ifelse(strand(lnc_par) == "+", b, a)
weak_over_up <- weak_overlap & over_up
weak_over_down <- weak_overlap & over_down
# extract weakly overlapping lncRNA-mRNA pairs 
weakmRNAup <- mRNA_par[weak_over_up]
weakmRNAup$gene_id
weaklncup <- lnc_par[weak_over_up]
weaklncup$gene_id
weakmRNAdown <- mRNA_par[weak_over_down]
weakmRNAdown$gene_id
weaklncdown <- lnc_par[weak_over_down]
weaklncdown$gene_id
upID <- cbind(weaklncup$gene_id,weakmRNAup$gene_id)
colnames(upID) <- c("up_lncRNA","mRNA")
downID <- cbind(weaklncdown$gene_id,weakmRNAdown$gene_id)
colnames(downID) <- c("down_lncRNA","mRNA")
write.table(downID,"downlncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")
write.table(upID,"uplncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")

# Switch the strand of transcripts:
lnc_sw <- lnc_rna
strand(lnc_sw) <- ifelse(strand(lnc_rna) == "+", "-", "+")
lnc_sw_start <- resize(lnc_sw, 1, "end")
# Detect DNC transcripts (start within 500 bp from TSS of a known gene):
dnc_intervals <- suppressWarnings(trim(flank(mRNA, 500)))
# find divergent antisense lncRNA and mRNA pairs
DNC_hits <- findOverlaps(lnc_sw_start,dnc_intervals)
lnc_DNC <-lnc_sw_start[queryHits(DNC_hits)]
lnc_DNC$gene_id
mRNA_DNC <- dnc_intervals[subjectHits(DNC_hits)]
mRNA_DNC$gene_id
divergentID <- cbind(lnc_DNC$gene_id,mRNA_DNC$gene_id)
colnames(divergentID) <- c("divanti_lncRNA","mRNA")
write.table(divergentID,"divergent_antilncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")

# Find overlaps with genes in the antisense orientation:
over_as <- lnc_sw %over% mRNA
hits2 <- findOverlaps(lnc_sw, mRNA)
lnc_sw_start_par <- lnc_sw_start[queryHits(hits2)]
mRNA_par2 <- mRNA[subjectHits(hits2)]
# Calculate relative distance between sTSS and asTSS:
mRNA_par2_start <- resize(mRNA_par2, 1, "start")
dist <- width(pgap(lnc_sw_start_par, mRNA_par2_start))
rel_dist <- dist  / width(mRNA_par2)
# Find convergent, pas-AS and distal AS transcripts:
conv <- rel_dist <= 0.5
pas_as <- rel_dist > 0.5 & rel_dist <= 1.2
distal <- rel_dist > 1.2
convmRNA <- mRNA_par2[conv]
convmRNA$gene_id
conlncRNA <- lnc_sw_start_par[conv]
conlncRNA$gene_id
pasAs_mRNA <- mRNA_par2[pas_as]
pasAs_mRNA$gene_id
pasAs_lncRNA <- lnc_sw_start_par[pas_as]
pasAs_lncRNA$gene_id
DistalmRNA <- mRNA_par2[distal]
DistalmRNA$gene_id
Distal_asRNA <- lnc_sw_start_par[distal]
Distal_asRNA$gene_id

convergentID <- cbind(conlncRNA$gene_id,convmRNA$gene_id)
colnames(convergentID) <- c("convanti_lncRNA","mRNA")
pasID <- cbind(pasAs_lncRNA$gene_id,pasAs_mRNA$gene_id)
colnames(pasID) <- c("pasanti_lncRNA","mRNA")
distalID <- cbind(Distal_asRNA$gene_id,DistalmRNA$gene_id)
colnames(distalID) <- c("distalanti_lncRNA","mRNA")
write.table(convergentID,"convergent_antilncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")
write.table(pasID,"pas_antilncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")
write.table(distalID,"distal_antilncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")
totaln_mRNA <- rbind(genicID,upID,downID,convergentID,pasID,distalID,divergentID)
#except intergenic lncRNA-mRNA pairs, we keep others
write.table(totaln_mRNA,"lncRNA_mRNA.txt",quote =F,row.names = F,sep = ",")
#generate interlncRNA_mRNA pairs
#generate interlncRNA Granges object
# scripts from sebastian group to annotate noncoding transcripts
annotateNoncodingTranscripts  <- function(tx, genes) {
  library(GenomicRanges)
  tx <- granges(tx)
  genes <- granges(genes)
  if (!identical(seqinfo(tx), seqinfo(genes))) { stop("Ensure identical seqinfo!") }
  # Prepare the output data frame:
  mat <- matrix(data = FALSE, nrow = length(tx), ncol = 8, dimnames = list(NULL, c("Genic", "Up", "Down", "Divergent", "Convergent", "TTS_AS", "Distal_AS", "Intergenic")))
  # Find strong overlaps on the sense strand (either transcript or gene is covered by at least 50% width):
  over_s <- tx %over% genes
  hits <- findOverlaps(tx, genes)
  tx_par <- tx[queryHits(hits)]
  genes_par <- genes[subjectHits(hits)]
  overlap <- width(pintersect(tx_par, genes_par))
  over_tx <- overlap / width(tx_par)
  over_genes <- overlap / width(genes_par)
  strong_overlap <- over_tx >= 0.5 | over_genes >= 0.5
  mat[over_s, "Genic"] <- as.logical(unlist(tapply(strong_overlap, queryHits(hits), any)))
  # Detect weak upstream and downstream overlaps on the same strand:
  weak_overlap <- !strong_overlap
  a <- start(tx_par) <= start(genes_par) & end(tx_par) < end(genes_par)
  b <- start(tx_par) > start(genes_par) & end(tx_par) >= end(genes_par)
  over_up <- ifelse(strand(tx_par) == "+", a, b)
  over_down <- ifelse(strand(tx_par) == "+", b, a)
  weak_over_up <- weak_overlap & over_up
  weak_over_down <- weak_overlap & over_down
  mat[over_s, "Up"] <- as.logical(unlist(tapply(weak_over_up, queryHits(hits), any)))
  mat[over_s, "Down"] <- as.logical(unlist(tapply(weak_over_down, queryHits(hits), any)))
  # Switch the strand of transcripts:
  tx_sw <- tx
  strand(tx_sw) <- ifelse(strand(tx) == "+", "-", "+")
  tx_sw_start <- resize(tx_sw, 1, "end")
  # Detect DNC transcripts (start within 500 bp from TSS of a known gene):
  dnc_intervals <- suppressWarnings(trim(flank(genes, 500)))
  mat[, "Divergent"] <- tx_sw_start %over% dnc_intervals
  # Find overlaps with genes in the antisense orientation:
  over_as <- tx_sw %over% genes
  hits2 <- findOverlaps(tx_sw, genes)
  tx_sw_start_par <- tx_sw_start[queryHits(hits2)]
  genes_par2 <- genes[subjectHits(hits2)]
  # Calculate relative distance between sTSS and asTSS:
  genes_par2_start <- resize(genes_par2, 1, "start")
  dist <- width(pgap(tx_sw_start_par, genes_par2_start))
  rel_dist <- dist  / width(genes_par2)
  # Find convergent, TTS-AS and distal AS transcripts:
  conv <- rel_dist <= 0.5
  tts_as <- rel_dist > 0.5 & rel_dist <= 1.2
  distal <- rel_dist > 1.2
  mat[over_as, "Convergent"] <- as.logical(unlist(tapply(conv, queryHits(hits2), any)))
  mat[over_as, "TTS_AS"] <- as.logical(unlist(tapply(tts_as, queryHits(hits2), any)))
  mat[over_as, "Distal_AS"] <- as.logical(unlist(tapply(distal, queryHits(hits2), any)))
  # Transcripts which are FALSE in all these classes are intergenic:
  mat[, "Intergenic"] <- ifelse(rowSums(mat[, 1:7]) > 0, FALSE, TRUE)
  return(mat)
}
# import adjusted araport11 gff file from sebastian group 
genes_araport <- read.csv("genes_araport_adj.bed",sep="\t",quote = "",header = F, stringsAsFactors = F,dec=".")
colnames(genes_araport) <- c("seqnames","start","end","gene_id","", "strand","thickstart","thickend")
# open adjusted araport11 Large Granges from desktop
genes_araport_adj <- readRDS("genes_araport_adj.RDS")
# extract lncRNAs/antilncRNA/mRNA and their IDs from  genes_araport_adj
senslncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="lnc_RNA"),]
antilncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="antisense_lncRNA"),]
lnc_rna <- c(senslncRNA,antilncRNA)
lncRNAID <- c(senslncRNA$gene_id,antilncRNA$gene_id)
lncRNA<- genes_araport[genes_araport$gene_id %in% lncRNAID, ]
gr_mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]

# run the scripts from sebastian group
annotateNoncodingTranscripts(gr_lncRNA,gr_mRNA)
#define the transcripts result
classlncRNA <- annotateNoncodingTranscripts(gr_lncRNA,gr_mRNA)
#view it and make it as a data.frame format
classlncRNA <- annotateNoncodingTranscripts(gr_lncRNA,gr_mRNA)
View(classlncRNA)
classlncRNA <- as.data.frame(classlncRNA) 
#sort  lncRNA
inter_lncRNA <-lncRNA[which(classlncRNA$Intergenic),]
## generate lincRNA_mRNA pairs
# others are intergenic non-coding RNAs and mRNAs
# load mRNA and intergenic lncRNA genome information
linter <- read.table("interlncRNA.txt", sep="\t", quote="", header=F, stringsAsFactors=F)
linter$V3 <- linter$V3+1000
MRNA <- read.table("mRNA.txt", sep="\t", quote="", header=F, stringsAsFactors=F)
# flip the strand of lncRNA in order to look for its antisense strand of codingRNA
linter$V6 <- ifelse(linter$V6 == "+", "-", "+")
linter$genes <- rep(NA,dim(linter)[1])
#A loop for production of lincRNA-Mrna PAIRS
for (i in 1:dim(linter)[1]) {
  genename<-c()
  for (j in 1:dim(MRNA)[1]) {
    # if we look at the same gene:
    if (MRNA$V1[j]==linter$V1[i]) {
      if (MRNA$V6[j]==linter$V6[i]) {
        #create the sequences for interval i and gene j
        X<-seq(linter$V2[i],linter$V3[i])
        Y<-seq(MRNA$V2[j],MRNA$V3[j])
        # if interval and gene colocalize = if any of Y belongs to X:
        if(any(Y %in% X)){
          # add gene j name to genename
          genename<-paste0(MRNA$V4[j]," ",genename)
        }
        #if there is gene names in genename, add this names to the genes column of $
        if(length(genename)>0){linter$genes[i]<-genename}
      }
    }
    
  }
}

lincRNA_mRNA <- cbind(linter$V4,linter$genes)
lincRNA_mRNA <- na.omit(lincRNA_mRNA)
lincRNA_mRNA[180,] <- c("AT2G09160","AT2G40710")
lincRNA_mRNA[279,] <- c("AT3G08985","AT3G57785")
lincRNA_mRNA <- as.data.frame(lincRNA_mRNA)
## remember there are some mistakes in the list. check carefully and save the table

write.table(lincRNA_mRNA,file = "lincRNA_mRNA",quote = F,sep=" ",dec=".",row.names = T,col.names = T)
# because the format of lincRNA_mRNA generated is weird, I reorganized in the text file and import it into R
lincRNA_mRNA <- read.table("lincRNA_mRNA.txt", sep="\t", quote="", header=F, stringsAsFactors=F)
