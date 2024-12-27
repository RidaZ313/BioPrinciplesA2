# Note: Used code from course materials. 

# Directory path
path = "C:\\Users\\ridaz\\OneDrive\\Desktop\\UCD\\Autumn 2024\\ANAT40040 - Biological Prinicples and Cellular Organisation\\Assignment 2"

# untar folder.
folder_name = "brca_tcga_pan_can_atlas_2018.tar.gz"
folder = paste(path, folder_name, sep = "/")
untar(folder)

# Go to new path
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)

# Read files
data_mrna = read.delim("data_mrna_seq_v2_rsem.txt")
data_clinical = read.delim("data_clinical_patient.txt")
data_cna = read.delim("data_cna.txt")

# Removing column 1 and column 2
assay = as.matrix(data_mrna[,-c(1,2)])
rownames(assay) = data_mrna[,1]

# Build metadata.
metadata = matrix(0, dim(assay)[2],1)

# Find patient ids
pat_ids = data_clinical[,1]
pat_ids = pat_ids[-c(1:4)]

# Find the column in data_cna containing the CNA level of ERBB2+
erbb2_row = which(data_cna$Hugo_Symbol == "ERBB2")

# Build metadata
for (i in 1:dim(assay)[2]){
  pat_barcode_initial = colnames(assay)[i]
  pat_barcode = substr(pat_barcode_initial, 1, 12)
  pat_barcode = gsub("\\.", "-",pat_barcode)
  
  idx = which(pat_barcode == pat_ids)
  idx2 = grep(pat_barcode_initial, colnames(data_cna))
  
  if (length(idx2) != 0){
    metadata[i,1] = 1*(as.numeric(data_cna[erbb2_row, idx2])>0)
  }
}

metadata[is.na(metadata)] = 0

colnames(metadata) = c("ERBB2")

#Install biocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2
BiocManager::install("DESeq2")

library(DESeq2)

# Build DESeq Object
assay[is.na(assay)] = 0  # Impute with zeros the NA
assay[assay<0] = 0

# Filter data
smallestGroupSize = 3
keep = rowSums(assay >= 10) >= smallestGroupSize
assay = assay[keep,]


dds <- DESeqDataSetFromMatrix(countData = round(assay),
                              colData = metadata,
                              design = ~ ERBB2)

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

res = results(dds)

# print Top 10 most differentially expressed
res[order(res$padj)[1:10],]

vsd = vst(dds)
par(mfrow = c(1, 2))
plotPCA(vsd, intgroup=c("ERBB2"))

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("enrichplot", quietly = TRUE))
  install.packages("enrichplot")

# Add required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res_sig = res[res$padj<0.05,]

DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

head(DE_over)
head(DE_under)

go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(go_results_over))

go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(go_results_under))

# Create dotplot
dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")

if (!requireNamespace("pathview", quietly = TRUE))
  BiocManager::install("pathview")

if (!requireNamespace("ReactomePA", quietly = TRUE)) 
  BiocManager::install("ReactomePA", force = TRUE)

library(ReactomePA) # Might not be available in this version of R
library(pathview)


gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_results_over))

# Create dot plot
dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

print(head(kegg_results_under))

# Create dot plot
dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")

# # Error -  library ReactomePA might not be available in this version of R
# reactome_results_over =  enrichPathway(
#   gene          = gene_entrez_over[,2],
#   organism      = "human",   
#   pAdjustMethod = "BH",
#   pvalueCutoff  = 0.05,
#   qvalueCutoff  = 0.05,
# )
# 
# # Error -  library ReactomePA might not be available in this version of R
# reactome_results_under =  enrichPathway(
#   gene          = gene_entrez_under[,2],
#   organism      = "human",   
#   pAdjustMethod = "BH",
#   pvalueCutoff  = 0.05,
#   qvalueCutoff  = 0.05,
# )
# 
# 
# print(head(reactome_results_over))
# 
# # Create dot plot
# dotplot(reactome_results_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")
# 
# print(head(reactome_results_under))
# 
# # Create dotplot
# dotplot(reactome_results_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")

go_results_under_pw = pairwise_termsim(go_results_under)
treeplot(go_results_under_pw)+ ggtitle("GO Enrichment Under Expressed")

kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")

top_DE = order(res$padj)

vsd_DE = assay(vsd)[top_DE[1:20],]

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

# Add package
library(pheatmap)

annotation_colors = list(ERBB2 = c(ERBB2 = "#1f78b4", Less = "#33a02c"))

annotation_col = data.frame(ERBB2 = as.matrix(metadata[,1]))
rownames(annotation_col) = colnames(vsd)

# Create heatmap
pheatmap(
  vsd_DE,
  cluster_rows = TRUE,      
  cluster_cols = TRUE,  
  scale = 'row',
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col)



