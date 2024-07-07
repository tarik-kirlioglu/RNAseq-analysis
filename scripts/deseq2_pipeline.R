#load to libraries

library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(tidyverse)

#load to counts data
counts <- read.csv("counts_celldensity.csv", header = T, row.names = 1)

#create to col data
col_data <- data.frame(condition=c("high", "high", "low", "low"))
row.names(col_data) <- colnames(counts)

#check rownames equal to colnames
all(colnames(counts) == rownames(col_data))

#construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~condition)

#filtering to low counts
kepp <- rowSums(counts(dds)) >=10
dds <- dds[kepp,]

#run DESeq2
dds <- DESeq(dds)

#retrival to results
res <- results(dds, contrast = c("condition", "high", "low"))

#retrival to normalized counts
normalized_counts <- counts(dds, normalized = T)

#retrival gene symbol

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

#pca plot
vstdata <- vst(dds,
               blind = F)
plotPCA(vstdata,
        intgroup="condition")

#plotMA
plotMA(res, alpha=0.05)

#plotMA with ggplot
res$significant <- ifelse(res$pvalue <= 0.05, "yes", "no")
ggplot(res, aes(x=log(baseMean), y=log2FoldChange, color= significant))+
  geom_point()

#volcano plot
EnhancedVolcano(res, x="log2FoldChange", y="padj", lab = res$symbol)

#heatmap with significant normalized counts
signi <- as.data.frame(res) %>% 
  subset(subset=res$pvalue <= 0.05)

all_signi <- merge(normalized_counts, signi, by=0)
counts_signi <- all_signi[,2:5]
row.names(counts_signi) <- all_signi$Row.names

pheatmap(log2(counts_signi + 1), 
         scale = "row", 
         show_rownames = F, 
         treeheight_row = 0,
         treeheight_col = 0)

#heatmap with top 20 genes 
top20 <- head(all_signi[order(all_signi$pvalue),],20)
top20_counts <- top20[ ,2:5]
row.names(top20_counts) <- top20$symbol
head(top20_counts)
pheatmap(log2(top20_counts + 1))

#save the results
write.csv(res, "deseq_cell_density.csv")
write.csv(normalized_counts, "normalize_counts.csv")
