#STEP 1: load to libraries

library(DESeq2)
library(tidyverse)
library(ggplot2)

#STEP 2: load to datasets

counts_data <- read.csv("counts.csv",
                        header = T,
                        sep =",",
                        row.names = 1)

col_data <- data.frame(condition=c("control", "control", "treated", "treated"))
rownames(col_data) <- colnames(counts)

#check rownames equal colnames

all(colnames(counts_data)==rownames(coldata))

# STEP 3: construct DESeqDataSet

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design =  ~condition)

#STEP 4: filtering to low counts and analyze PCA 

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

vstdata <- vst(dds,
               blind = F)
plotPCA(vstdata,
        intgroup="condition")

#STEP 5: set to level

dds$condition <- relevel(dds$condition,
                         ref = "control")

#STEP 6: run DESeq, plotMA and write to data

dds <- DESeq(dds)
plotMA(dds)
res <- results(dds, tidy = T)
write.csv(res, "log_fold_change.csv")
