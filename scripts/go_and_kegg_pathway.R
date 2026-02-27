#importing libraries
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(org.Hs.eg.db)

#reading results
res <- read.csv("res_and_normalized.csv", header = T, row.names = 2)

#selecting significant genes
significant_genes <- res %>% 
  filter(padj <= 0.05 & abs(log2FoldChange) > 1)

#getting ensembl ids
genelist <- rownames(significant_genes)

#enrichment gene ontology analysis, ontology(ont) can be changed, "MF", "BP, and "CC" can be entered.
ego <- enrichGO(gene = genelist,
                OrgDb = "org.Hs.eg.db",
                keyType = "ENSEMBL",
                ont = "BP",
                readable = T)

#plotting ego
barplot(ego)

#getting uniprot ids with fold change values for pathview
gene_data <- significant_genes$log2FoldChange
names(gene_data) <- significant_genes$uniprot_id
uniprotid <- significant_genes$uniprot_id

#searching kegg organism code
homo_sapiens <- search_kegg_organism("Homo sapiens", by = "scientific_name")

#enrichment kegg pathway analysis
ekegg <- enrichKEGG(gene = uniprotid,
                   keyType = "uniprot",
                   organism = "hsa")

#visualization
browseKEGG(ekegg, pathID = "hsa04080")

hsa04080 <- pathview(gene.data = gene_data,
                     pathway.id = "hsa04080",
                     species = "hsa",
                     gene.idtype = "uniprot")
