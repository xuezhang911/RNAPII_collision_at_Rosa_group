## using downloaded feature counts to continue gene expression differentiation analysis
# data is from GSE135607, 1week whole seedling, the trimmed reads were aligned to the arabidopsis genome Tair10 using topha2
#50um ABA treatment for 3 hours... some ABA markers increase, CDPK10 slightly increase but NRT2.9 doesn't change in WT. 
data <- read.table("GSE135607_JKZ218_featureCounts_output.txt",header = T)
# only use ABA treatment data in col0
sampleNames <- c("WTc_r1","WTc_r2","WTc_r3","WTABA_r1","WTABA_r2","WTABA_r3")
# we only need count data 
 names(data)[7:12] <- sampleNames
 countData <- as.matrix(data[7:12])
 rownames(countData) <- data$Geneid
 database <- data.frame(name=sampleNames, condition=c("control","control","control","ABA","ABA","ABA"))
 rownames(database) <- sampleNames 
 # load library deseq2
 library(DESeq2)
 ## generate dds can start from summarizedexperiment object or a (count matrix and a sample information table by featurecount)
 dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
 dds <- dds[ rowSums(counts(dds)) > 1, ]
 # generate DGElist
 dds <- DESeq(dds)
 res <- results(dds) # diffirential gene expression matrix
 # or same res <- DESeq2::results(dds)
 # transform to dataframe
 df <- as.data.frame(res)
 id <- row.names(df)
 df <- cbind(id,df)
 df <- df[order(df$log2FoldChange,decreasing = T),]
 ## including single experiment data 
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
 resdata <- resdata[order(resdata$log2FoldChange,decreasing = T),]
# vocalno plot
 library(ggplot2)
 resdata$change <- as.factor(
   ifelse(
     resdata$padj<0.01 & abs(resdata$log2FoldChange)>1,
     ifelse(resdata$log2FoldChange>1, "Up", "Down"),
     "NoDiff"
   )
 )
 valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
   geom_point(alpha=0.8, size=1) + 
   theme_bw(base_size=15) + 
   theme(
     panel.grid.minor=element_blank(),
     panel.grid.major=element_blank()
   ) + 
   ggtitle("DESeq2 Valcano") + 
   scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
   geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
   geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
 valcano
 # PCA plot
 rld <- rlog(dds)
 pcaData <- plotPCA(rld, intgroup=c("condition", "name"), returnData=T)
 percentVar <- round(100*attr(pcaData, "percentVar"))
 pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=name)) + 
   geom_point(size=3) + 
   ggtitle("DESeq2 PCA") + 
   xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
   ylab(paste0("PC2: ", percentVar[2], "% variance"))
 pca
 
 library(pheatmap)
 
 select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
 nt <- normTransform(dds)
 log2.norm.counts <- assay(nt)[select,]
 df <- as.data.frame(colData(dds)[, c("name", "condition")])
 pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
 
# use the variance stablizing transcformation implemented with the vst function 
vsd <- DESeq2::vst(dds)
class(vsd)
head(colData(vsd), 3)
# can also proce PCA Plot # not super beautiful
DESeq2::plotPCA(vsd, intgroup=c("condition", "name"))
# MDS plot from edgR package
hist(res$pvalue)
res$log10BaseMean <- log10(res$baseMean)
res$mlog10PValue <- -log10(res$pvalue)
rowData(dds)$log10Dispersion <- log10(rowData(dds)$dispersion)
rowData(dds)$DESeq2_dex_trt_vs_untrt <- res
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold = 1)
summary(resLFC1)
table(resLFC1$padj < 0.1)
plotCounts(dds, gene = "AT1G18880", intgroup = "condition", 
           normalized = TRUE, transform = FALSE)
plotCounts(dds, gene = "AT1G18890", intgroup = "condition", 
           normalized = TRUE, transform = FALSE)
# MA plot: useful overview for an experiment with two-group comparision 
DESeq2::plotMA(res, ylim = c(-5, 5))
# visualize the most significant genes  very useful
suppressPackageStartupMessages({
  library(pheatmap)
})
mat <- assay(vsd)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
dfx <- as.data.frame(colData(vsd)[, c("condition","name")])
pheatmap(mat, annotation_col = dfx)
# exporting gene anotataed results to csv file
library(org.At.tair.db)
res$symbol <- mapIds(org.At.tair.db, keys = gsub("\\.[0-9]+$", "", row.names(res)),
       column = c('SYMBOL'), keytype = 'TAIR',multiVals = "first")
resOrdered <- res[order(res$padj), ]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[seq_len(100), ]
write.table(cbind(id = rownames(resOrderedDF), resOrderedDF), 
            file = "results.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)