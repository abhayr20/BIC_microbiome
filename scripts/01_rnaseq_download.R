# Load packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

path = '/home/user/Desktop/abhay/BIC_microbiome/'
setwd(path)

# get a list of projects
GDCprojects = getGDCprojects()
getProjectSummary('TCGA-BRCA')

# build a query to retrieve gene expression data ------------
query_brca <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')

query_clinical <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML")

getResults(query_brca)
getResults(query_clinical)

# download data - GDCdownload
GDCdownload(query_brca)

# prepare data
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
saveRDS(object = tcga_brca_data, file = "tcga_brca.RDS", compress = FALSE)
tcga_brca_data = readRDS(file = "tcga_brca.RDS")
