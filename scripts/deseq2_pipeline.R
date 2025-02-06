#importing libraries
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pheatmap)

#load to counts data
counts <- read.csv("gene_counts.csv", header = T, row.names = 1)

#create to col data
colData <- data.frame(condition=c("control", "control", "control", "knockout", "knockout", "knockout"))
rownames(colData) <- colnames(counts)

#check rownames equal to colnames
all(rownames(colData) == colnames(counts))

#construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~condition)

#filtering to low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#setting reference level
dds$condition <- relevel(dds$condition, ref = "control")

#run DESeq2
dds <- DESeq(dds)

#retrival to results
res <- results(dds)

#retrival to normalizedcounts
normalized_counts <- counts(dds, normalized = T)

#retrival uniprot id and  gene symbol
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

res$uniprot_id <- mapIds(org.Hs.eg.db,
                         keys = rownames(res),
                         keytype = "ENSEMBL",
                         column = "UNIPROT")

res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

#convert res to dataframe and merge res and normalized counts
res <- as.data.frame(res)
res_and_normalized <- merge(res, normalized_counts, by=0)

#plotMA graph
plotMA <- plotMA(dds, alpha=0.05)

#volcano plot
enhancedvolcano <- EnhancedVolcano(res_and_normalized,
                                   x="log2FoldChange",
                                   y="padj",
                                   lab = res_and_normalized$symbol)

#retrieving significant normalized counts
signi_normalized_counts <- res_and_normalized %>% 
  filter(res$padj <= 0.05 & abs(res$log2FoldChange) > 1) %>% 
  dplyr::select(10:15)

#heatmap with significant normalized counts
heatmap_normalized <- pheatmap(log2(signi_normalized_counts + 1), 
         scale = "row", 
         show_rownames = F, 
         treeheight_row = 0,
         treeheight_col = 0)

#retrieving significant top20 genes
top20_genes <- res_and_normalized %>% 
  filter(res$padj <= 0.05 & abs(res$log2FoldChange) > 1) %>% 
  arrange(padj) %>% 
  dplyr::select(9:15) %>% 
  head(20)
  
rownames(top20_genes) <- top20_genes$symbol
top20_genes$symbol <- NULL

#heatmap with top20 genes
heatmap_top20genes <- pheatmap(log2(top20_genes + 1), 
         scale = "row")

#save the results
write.csv(res, "deseq.csv")
write.csv(normalized_counts, "normalized_counts.csv")
write.csv(res_and_normalized, "res_and_normalized.csv")
