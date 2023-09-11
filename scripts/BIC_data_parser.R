library(tidyverse)
library(SummarizedExperiment)
library(GenomicDataCommons)
library(TCGAbiolinks)
library(TCGAutils)
library(magrittr)

ls("package:TCGAutils")
setwd("/data/abhay/brca")

###################################################################
# Helper functions
all_unique <- function(input_vector) {
  length(input_vector) == length(unique(input_vector))
}

all_unique_vector <- function(input_vector) {
  return(length(input_vector) == length(unique(input_vector)))
}

truncate_strings <- function(input_vector, num_characters = 12) {
  truncated_vector <- substr(input_vector, 1, num_characters)
  return(truncated_vector)
}

# To parse the phylogenetic info from metagenomic cnt_mtx
extract_parts <- function(strings) {
  result <- data.frame(otu = sprintf("Otu%04d", 1:nrow(brca_kraken)), 
                       kingdom = character(length(strings)),
                       phylum = character(length(strings)),
                       class = character(length(strings)),
                       order = character(length(strings)),
                       family = character(length(strings)),
                       genus = character(length(strings)))
  
  for (i in seq_along(strings)) {
    parts <- unlist(strsplit(strings[i], "\\."))  # Split the string by periods
    
    for (part in parts) {
      if (grepl("^k__", part)) result[i, "kingdom"] <- gsub("k__", "", part)
      if (grepl("^p__", part)) result[i, "phylum"] <- gsub("p__", "", part)
      if (grepl("^c__", part)) result[i, "class"] <- gsub("c__", "", part)
      if (grepl("^o__", part)) result[i, "order"] <- gsub("o__", "", part)
      if (grepl("^f__", part)) result[i, "family"] <- gsub("f__", "", part)
      if (grepl("^g__", part)) result[i, "genus"] <- gsub("g__", "", part)
    }
  }
  
  return(result)
}
###################################################################

###################################################################
# Read clinical, metagenomic and transcriptomic data
tcga_clinical <- read_tsv("data/TCGA-BIC-clinical.tsv")                      # clinical parameters 
tcga_meta     <- read.csv("data/Metadata-TCGA-Kraken-17625-Samples.csv")     # metagenomic metadata, X is a unique patient sample
tcga_kraken   <- read.csv("data/Kraken-TCGA-Raw-Data-17625-Samples.csv")     # metagenomic cnt_mxt
tcga_rnaseq   <- readRDS(file = "scripts/tcga_brca.RDS")                     # BRCA RNAseq data object
###################################################################


###################################################################
# Filter clinical, metadata and metagenomic cnt_mtx to only BRCA
# Find patient IDs using TCGAutils
brca_meta <- tcga_meta %>%
  dplyr::filter(toupper(case_uuid) %in% toupper(tcga_clinical$`Other Patient ID`) &
           sample_type %in% c("Primary Tumor", "Metastatic") &
           investigation == "TCGA-BRCA" &
           experimental_strategy == "RNA-Seq")

brca_meta <- brca_meta[!duplicated(brca_meta$case_uuid), ] # Remove duplicate patient samples
brca_map <- UUIDtoBarcode(unique(brca_meta$case_uuid))     # Use TCGAutils to get patient IDs using case UUIDs
brca_clinical <- tcga_clinical[tcga_clinical$`Patient ID` %in% brca_map$submitter_id , ]

brca_map <- brca_map %>%
  left_join(brca_clinical[, c("Patient ID", "Diagnosis Age", "Race Category", "Sex", "Subtype")], by = c("submitter_id" = "Patient ID"))

colnames(brca_map) <- c('case_uuid', 'patient_id', 'age', 'race', 'sex', 'subtype')
brca_map$X <- brca_meta$X

brca_kraken <- tcga_kraken %>%
  dplyr::filter(X %in% brca_meta$X) %>%
  column_to_rownames("X") %>%
  t() %>%
  as.data.frame()
###################################################################


###################################################################
# Create taxonomic classification matrix by parsing metagenomic cnt_mtx rownames
microbes <- rownames(brca_kraken)
rownames(brca_kraken) <- sprintf("Otu%04d", 1:nrow(brca_kraken))
tax_mat <- extract_parts(microbes) 
###################################################################

###################################################################
# RNAseq cnt_mtx
'''
Available assays in SummarizedExperiment : 
=> unstranded
=> stranded_first
=> stranded_second
=> tpm_unstrand      ***
=> fpkm_unstrand
=> fpkm_uq_unstrand
'''

# Read rnaseq cnt_mtx as a df
brca_rnaseq <- as.data.frame(assay(tcga_rnaseq, 'tpm_unstrand')) 

# Filter rows (reads): Change ENSG to HGSV, remove duplicate HGSV terms and zero reads
map_genes <- data.frame(ENGS=brca_data@rowRanges@elementMetadata@listData$gene_id,
                       HGSV=brca_data@rowRanges@elementMetadata@listData$gene_name)

unq_genes <- map_genes[!(duplicated(map_genes$HGSV)),]
brca_rnaseq <- brca_rnaseq[(rownames(brca_rnaseq) %in% unq_genes$ENGS),]
rownames(brca_rnaseq) <- unq_genes$HGSV

zero_reads <- brca_rnaseq[rowSums(brca_rnaseq[]) == 0,]
brca_rnaseq <- brca_rnaseq[rowSums(brca_rnaseq[]) > 0,]

#Filter cols (patients): Remove non-tumor samples and duplicate
map_pts <- data.frame(
  patient_id = tcga_rnaseq@colData@listData$patient,
  barcode = tcga_rnaseq@colData@listData$barcode,
  sample_type = tcga_rnaseq@colData@listData$definition
)

map_pts <- map_pts %>%
  dplyr::filter(sample_type %in% c("Metastatic", "Primary solid Tumor")) %>%
  dplyr::filter(patient_id %in% brca_map$patient_id)

# 23 patients where removed which have more than 1 sample
dup_pts <- map_pts[duplicated(map_pts$patient_id),] 
unq_pts <- map_pts[!duplicated(map_pts$patient_id),]

brca_map <- merge(brca_map, unq_pts, by = 'patient_id') 
brca_rnaseq <- brca_rnaseq[,unq_pts$barcode]            

colnames(brca_rnaseq) <- unq_pts$patient_id # 
index <- match(colnames(brca_rnaseq), brca_map$patient_id)
brca_map <- brca_map[index,]
colnames(brca_rnaseq) <- brca_map$X

###################################################################

######################################################
# Check and compare case UUIDs
v1 <- toupper(brca_clinical$`Other Patient ID`)
v2 <- toupper(brca_map$case_id)
matching_elements <- intersect(v1, v2)
unique_to_v1 <- setdiff(v1, v2)
unique_to_v2 <- setdiff(v2, v1)
num_duplicated_v1 <- sum(duplicated(v1))
num_duplicated_v2 <- sum(duplicated(v2))
duplicated_v2 <- v2[duplicated(v2)]

#Final check
alpha <- read.csv("/data/abhay/brca/results/brca_alpha_div.csv", row.names = 1)
all(sort(colnames(brca_rnaseq)) == sort(rownames(alpha)))
######################################################

# Write parsed files
write.csv(brca_kraken, file = "cleaned_data/otu_mat.csv", row.names = TRUE)
write.csv(brca_meta, file = "cleaned_data/sample_mat.csv", row.names = FALSE)
write.csv(tax_mat, file = "cleaned_data/BIC_tax_mat.csv", row.names = FALSE)
write.csv(brca_map, file = "cleaned_data/brca_map.csv", row.names = FALSE)
write.csv(brca_clinical, file = "cleaned_data/brca_clinical.csv", row.names = FALSE)
write.csv(brca_rnaseq, file = "cleaned_data/brca_rnaseq.csv" , row.names = TRUE)

