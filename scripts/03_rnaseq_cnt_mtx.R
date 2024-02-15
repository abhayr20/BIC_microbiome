# Load packages
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(readr)

path = '/home/user/Desktop/abhay/BIC_microbiome/'
setwd(path)

brca_data = readRDS(file = "tcga_brca.RDS")
m <- assay(brca_data, 'tpm_unstrand')
m <- as.data.frame(m)

# Change ENSG to HGSV, remove duplicate HGSV terms and zero reads
gene_map <- data.frame(ENGS=brca_data@rowRanges@elementMetadata@listData$gene_id,
                       HGSV=brca_data@rowRanges@elementMetadata@listData$gene_name)

dup_map <- gene_map[duplicated(gene_map$HGSV),]
dup_m <- m[rownames(m) %in% dup_map$ENGS,]

unq_map <- gene_map[!(duplicated(gene_map$HGSV)),]
m <- m[(rownames(m) %in% unq_map$ENGS),]
rownames(m) <- unq_map$HGSV

zero_reads <- m[rowSums(m[]) == 0,]
m <- m[rowSums(m[]) > 0,]

# Save final count matrix
write.csv(m, file = "cleaned_data/brca_cnt_mtx.csv", row.names = TRUE)
