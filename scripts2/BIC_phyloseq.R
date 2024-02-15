library("tidyverse")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble") 
library("phyloseq")

setwd("/data/abhay/brca")

######################################################
# Phyloseq installation and documentation
# library(BiocManager)
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

bic <- phyloseq(OTU, TAX, samples)


# Normalize number of reads in each sample using median sequencing depth.
# The number of reads used for normalization is **`r sprintf("%.0f", total)`**. 
total = median(sample_sums(bic))
standf = function(x, t=total) round(t * (x / sum(x)))
bic = transform_sample_counts(bic, standf)

# beta-diversity Clusterinig
ndms_distance <- ordinate(physeq=bic, method="NMDS", distance = "bray") #nmds
ordinate.df <- as.data.frame(ndms_distance$points)
ordinate.df$race <- as.factor(df_mat$race)

ggplot(ordinate.df) +
  geom_point(aes(x=MDS1, y=MDS2, col=race), size = 2)

pcoa_dist <- ordinate(physeq=bic, method="PCoA", distance = "bray") 
pcoa.df <- as.data.frame(pcoa_dist$vectors)
pcoa.df$race <- as.factor(df_mat$race)

beta_div <- ggplot(pcoa.df) +
  geom_point(aes(x=Axis.1, y=Axis.2, col=race), size = 2)
ggsave("results/beta_div.png", plot = beta_div)

# Alpha-diversity clustering
alpha_diversity<- estimate_richness(bic, split = TRUE, measures = NULL)
write.csv(alpha_diversity, file = "results/brca_alpha_div.csv", col.names = T, row.names = T)

# save.image(file = "rdata/brca_phyloseq.Rdata")
load("rdata/brca_phyloseq.Rdata")

brca_map2 <- read.csv("cleaned_data/brca_map.csv")
rownames(brca_map2) <- brca_map2$X

enrich <- as.data.frame(bic@otu_table) %>%
  dplyr::filter(rowMeans(.) >= 20) %>%
  t() %>%
  as.data.frame()  

enrich <- cbind(enrich, brca_map2[,4])
colnames(enrich)[246] <- "race"
enrich <- enrich[order(enrich$race),]

enrich_long <- pivot_longer(enrich, cols = c(1:245)) %>% View()
  group_by(race, name) %>% 
  summarise(rel_abundance = sum(value)) %>%
  View()

enrich_long <- enrich_long %>% arrange(desc(new_value)) 
enrich_long %>% group_by(race) %>% head(., n = 10) %>% View()

msleep %>%
  drop_na(sleep_rem, vore) %>%
  group_by(vore) %>%
  View()
