library(tidyverse)  
library(cluster)    
library(factoextra)
library(gridExtra) 

setwd("/Users/admin/Desktop/AIIMS/BIC_microbiome/")

tcga_kraken <- read.csv("data/Kraken-TCGA-Raw-Data-17625-Samples.csv")
tcga_meta <- read.csv("data/Metadata-TCGA-Kraken-17625-Samples.csv")
tcga_clinical_bic <- read_tsv("data/TCGA-BIC-clinical.tsv")

#Filtering method 1
tcga_meta_bic <- tcga_meta %>%
  filter(investigation == "TCGA-BRCA") %>%
  filter(sample_type == "Primary Tumor") %>%
  filter(experimental_strategy == "RNA-Seq")

#Filtering method 2
tcga_meta_bic2 <- tcga_meta_bic %>%
  filter(toupper(case_uuid) %in% toupper(tcga_clinical_bic$`Other Patient ID`))

# checking 
v1 <- toupper(tcga_clinical_bic$`Other Patient ID`)
v2 <- toupper(tcga_meta_bic2$case_uuid)
matching_elements <- intersect(v1, v2)
unique_to_v1 <- setdiff(v1, v2)
unique_to_v2 <- setdiff(v2, v1)
sum(duplicated(v1))
unique(duplicated(v2))
v2[duplicated(v2)]

# List of patients 
bic_list <- tcga_meta_bic2$X

tcga_kraken_bic <- tcga_kraken %>%
  filter(X %in% bic_list)

# Assuming you have a dataframe named 'my_dataframe'
write.csv(tcga_kraken_bic, file = "data/tcga_kraken_bic.csv", row.names = FALSE)
write.csv(tcga_kraken_bic, file = "data/tcga_meta_bic2.csv" , row.names = FALSE)

