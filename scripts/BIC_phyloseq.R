library(BiocManager)
library(phyloseq)

setwd("~/Desktop/AIIMS/BIC_microbiome")

######################################################
# Phyloseq installation and documentation
# BiocManager::install("phyloseq")
# browseVignettes("phyloseq") 
######################################################

# Load raw data
otu_mat <- read.csv("cleaned_data/otu_mat.csv", row.names = 1)      #Rows = species, cols = samples
tax_mat <- read.csv("cleaned_data/BIC_tax_mat.csv", row.names = 1)  #Rows = samples, cols = metadata
df_mat  <- read.csv("cleaned_data/sample_mat.csv", row.names = 1)   

# Transform to matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(df_mat)

carbom <- phyloseq(OTU, TAX, samples)
carbom
