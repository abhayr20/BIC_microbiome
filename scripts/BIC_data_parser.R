library(tidyverse)
setwd("/Users/admin/Desktop/AIIMS/BIC_microbiome/")

# Read data
tcga_kraken <- read.csv("data/Kraken-TCGA-Raw-Data-17625-Samples.csv")
tcga_meta <- read.csv("data/Metadata-TCGA-Kraken-17625-Samples.csv")
tcga_clinical_bic <- read_tsv("data/TCGA-BIC-clinical.tsv")

# Filter metadata with patient IDs
meta_bic <- tcga_meta %>%
  filter(toupper(case_uuid) %in% toupper(tcga_clinical_bic$`Other Patient ID`),
         sample_type %in% c("Primary Tumor","Metastatic"),
         investigation == "TCGA-BRCA",
         experimental_strategy == "RNA-Seq")

# Check and compare patient IDs
v1 <- toupper(tcga_clinical_bic$`Other Patient ID`)
v2 <- toupper(meta_bic$case_uuid)
matching_elements <- intersect(v1, v2)
unique_to_v1 <- setdiff(v1, v2)
unique_to_v2 <- setdiff(v2, v1)
num_duplicated_v1 <- sum(duplicated(v1))
num_duplicated_v2 <- sum(duplicated(v2))
duplicated_v2 <- v2[duplicated(v2)]

# Filter count matrix using BIC patient IDs
bic_list <- meta_bic$X
kraken_bic <- tcga_kraken %>%
  filter(X %in% bic_list)

# Organize the OTU matrix
kraken_bic <- kraken_bic %>% 
  column_to_rownames("X") %>%
  t() %>%
  as.data.frame()
microbes <- rownames(kraken_bic)
rownames(kraken_bic) <- sprintf("Otu%04d", 1:nrow(kraken_bic))

# Create taxonomic classification matrix
extract_parts <- function(strings) {
  result <- data.frame(otu = sprintf("Otu%04d", 1:nrow(kraken_bic)), 
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
tax_mat <- extract_parts(microbes)

# Write parsed files
write.csv(kraken_bic, file = "cleaned_data/otu_mat.csv", row.names = TRUE)
write.csv(meta_bic, file = "cleaned_data/sample_mat.csv", row.names = FALSE)
write.csv(tax_mat, file = "cleaned_data/BIC_tax_mat.csv", row.names = FALSE)
