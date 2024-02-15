library("tidyverse")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble") 
library("phyloseq")
library('vegan')
library('gridExtra')

path = '/home/user/Desktop/abhay/BIC_microbiome/'
setwd(path)


#Load EPIC immunescores
immune_scores <- read.csv("results/epic_immunedeconv_brca.csv", row.names = 1)
as.data.frame(t(immune_scores)) 
immune_scores <- immune_scores %>%
  rownames_to_column(var = 'samples') %>%
  select(samples, `T cell CD4+`, `T cell CD8+`, `NK cell`, `Macrophage`)
colnames(immune_scores) <- c('samples', 'CD4_T_cells', 'CD8_T_cells',
                             'NK_cells', 'Macrophages')


# Load OTU tables
#otu_mat <- read.csv("cleaned_data/otu_mat.csv", row.names = 1)      #Rows = species, cols = samples

#Removes 10 samples with less than 100k OTUs and 165 OTU with zero reads
otu_mat_filtered <- read.csv("cleaned_data/otu_mat.csv") %>%
  rename(otu = X) %>%
  pivot_longer(cols = -otu, names_to = 'samples') %>%
  pivot_wider(names_from = otu, values_from = value) %>%
  pivot_longer(-samples) %>%
  group_by(samples) %>%
  mutate(total = sum(value)) %>%
  filter(total > 100000) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(samples) %>%
  as.data.frame()

#Reshape otu_mat_filtered
rownames(otu_mat_filtered) <- otu_mat_filtered$samples
otu_mat_filtered <- otu_mat_filtered[, -1]
otu_mat <- as.data.frame(t(otu_mat_filtered))
rm(otu_mat_filtered)

#Load taxonomy data
tax_mat <- read.csv("cleaned_data/BIC_tax_mat.csv", row.names = 1)  #Rows = samples, cols = metadata

#Load Metadata
brca_map <- read_csv("cleaned_data/brca_map.csv") %>%
  select(X, subtype) %>%
  rename(samples = X)

metadata <- read_csv("cleaned_data/sample_mat.csv") %>%
  select(X, age_at_diagnosis, race) %>%
  rename(samples = X) %>%
  rename(age = age_at_diagnosis) %>%
  inner_join(., brca_map, by = 'samples') %>%
  mutate(subtype = str_replace_all(subtype, c("BRCA_Basal" = "Triple-neg",
                                              "BRCA_Normal" = "Triple-neg")),
         subtype = str_replace_all(subtype, "BRCA_LumA", "Luminal-A"),
         subtype = str_replace_all(subtype, "BRCA_LumB", "Luminal-B"),
         subtype = str_replace_all(subtype, "BRCA_Her2", "Her2-pos"),
         subtype = str_replace_na(subtype, "NA"),
         race    = str_replace_all(race, "WHITE", "Caucasian"),
         race    = str_replace_all(race, "ASIAN", "Asian"),
         race    = str_replace_all(race, "BLACK OR AFRICAN AMERICAN", "Black"),
         race    = str_replace_all(race, "AMERICAN INDIAN OR ALASKA NATIVE", "Native American"),
         race    = str_replace_all(race, "Not available", "NA")) 

rownames(metadata) <- metadata$samples
rm(brca_map)

# Transform to matrix
#otu_mat_filtered <- as.matrix(otu_mat_filtered)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
# samples = sample_data(metadata)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
bic <- phyloseq(OTU, TAX, metadata)

# Rarefy data
set.seed(19970911)
rarefy_even_depth(bic, sample.size = min(sample_sums(bic)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Normalize number of reads in each sample using median sequencing depth.
# total = median(sample_sums(bic))
# standf = function(x, t=total) round(t * (x / sum(x)))
# bic = transform_sample_counts(bic, standf)

# NMDS method for beta-diversity 
set.seed(1)
ndms_distance <- ordinate(physeq=bic, method="NMDS", distance = "bray")
ordinate.df <- as.data.frame(ndms_distance$points) %>%
  mutate(samples = rownames(.)) %>%
  inner_join(., metadata, by = 'samples') %>%
  filter(race != 'NA') %>%
  filter(race != "Native American") %>%
  filter(subtype != "NA")

her2  <- ordinate.df %>% filter(subtype == 'Her2-pos') 
tnbc  <- ordinate.df %>% filter(subtype == 'Triple-neg')
lum_a <- ordinate.df %>% filter(subtype == 'Luminal-A')
lum_b <- ordinate.df %>% filter(subtype == 'Luminal-B')
asian <- ordinate.df %>% filter(race == 'Asian')
black <- ordinate.df %>% filter(race == 'Black')
white <- ordinate.df %>% filter(race == 'Caucasian')

ggplot(ordinate.df) +
  geom_point(aes(x = MDS1, y = MDS2, col = race, shape = subtype), size = 2) +
  labs(title = "Beta-diversity clustering by Race and Cancer subtype")
ggsave("results/beta_diversity.png")

ggplot(ordinate.df) +
  geom_point(aes(x = MDS1, y = MDS2, col = race), size = 2) +
  labs(title = "Beta-diversity clustering by race") +
  facet_wrap(subtype~.)
ggsave("results/beta_diversity_race.png")


ggplot(ordinate.df) +
  geom_point(aes(x = MDS1, y = MDS2, col = subtype), size = 2) +
  labs(title = "Beta-diversity clustering by cancer subtype") +
  facet_wrap(race~.)
ggsave("results/beta_diversity_subtype.png")


# PCoA method for beta-diversity
pcoa_dist <- ordinate(physeq=bic, method="PCoA", distance = "bray") 
pcoa.df <- as.data.frame(pcoa_dist$vectors)%>%
  mutate(samples = rownames(.)) %>%
  inner_join(., metadata, by = 'samples')

beta_div <- ggplot(pcoa.df) +
  geom_point(aes(x=Axis.1, y=Axis.2, col=race), size = 2)
ggsave("results/beta_div.png", plot = beta_div)



# Alpha-diversity clustering
alpha_diversity<- estimate_richness(bic, split = TRUE, measures = NULL) 
# write.csv(alpha_diversity, file = "results/brca_alpha_div.csv", col.names = T, row.names = T)

alpha_immune <- alpha_diversity %>%
  rownames_to_column(var = 'samples') %>%
  select(samples, Shannon) %>%
  inner_join(., immune_scores, by = 'samples') %>%
  inner_join(., metadata, by = 'samples') %>%
  filter(race %in% c('Caucasian', "Black", "Asian")) %>%
  filter(subtype %in% c('Luminal-A', 'Luminal-B', 'Triple-neg', "Her2-pos"))

cor_shannon_CD4_T_cells <- cor.test(alpha_immune$Shannon, alpha_immune$CD4_T_cells, method = 'spearman')
cor_shannon_CD8_T_cells <- cor.test(alpha_immune$Shannon, alpha_immune$CD8_T_cells, method = 'spearman')
cor_shannon_macrophages <- cor.test(alpha_immune$Shannon, alpha_immune$Macrophages, method = 'spearman')

lm_shannon_CD4_T_cells_race <- lm(Shannon~CD4_T_cells + race, data=alpha_immune)
lm_shannon_CD8_T_cells_race <- lm(Shannon~CD8_T_cells + race, data=alpha_immune)
lm_shannon_macrophages_race <- lm(Shannon~Macrophages + race, data=alpha_immune)

lm_shannon_CD4_T_cells_subtype <- lm(Shannon~CD4_T_cells + subtype, data=alpha_immune)
lm_shannon_CD8_T_cells_subtype <- lm(Shannon~CD8_T_cells + subtype, data=alpha_immune)
lm_shannon_macrophages_subtype <- lm(Shannon~Macrophages + subtype, data=alpha_immune)

p <- paste("P-value:", round(cor_shannon_macrophages$p.value, digits=5))
rho <- paste("rho:", round(cor_shannon_macrophages$estimate, digits=2))
annotation <- paste(p, rho, sep="\n")

macrophage_race <- ggplot(alpha_immune, aes(x=Macrophages, y=Shannon, color=race))+
  geom_point() +
  geom_smooth(method="lm") +
  xlim(0, 0.20) +
  geom_text(aes(x=0.15, y=2, label=annotation), color="black", hjust = "left") +
  scale_color_manual(name=NULL,
                     values=c("blue", "black",  "red"),
                     breaks=c('Caucasian', "Black", "Asian"),
                     labels=c('Caucasian', "Black", "Asian")) +
  labs(title="There is no significant association between macropahge scores and Alpha diversity.\n\nThere is no significant difference between the three ethinicities.",
       x="Macrophage score calculated by EPIC",
       y="Shannon Diversity Index") +
  theme_classic()

macrophage_subtype <- ggplot(alpha_immune, aes(x=Macrophages, y=Shannon, color=subtype))+
  geom_point() +
  geom_smooth(method="lm") +
  xlim(0, 0.20) +
  geom_text(aes(x=0.15, y=2, label=annotation), color="black", hjust = "left") +
  scale_color_manual(name=NULL,
                     values=c("blue", "black",  "red", 'yellow'),
                     breaks=c('Luminal-A', 'Luminal-B', 'Triple-neg', "Her2-pos"),
                     labels=c('Luminal-A', 'Luminal-B', 'Triple-neg', "Her2-pos")) +
  labs(title="There is no significant difference between the four subtypes.",
       x="Macrophage score calculated by EPIC",
       y="Shannon Diversity Index") +
  theme_classic()

macrophage <- grid.arrange(macrophage_race, macrophage_subtype)
ggsave("results/macrophage_alpha_div.png", plot = macrophage, width = 8.5)

write.csv(alpha_diversity, file = "results/brca_alpha_div.csv", col.names = T, row.names = T)
#save.image(file = "rdata/brca_phyloseq.Rdata")
#load("rdata/brca_phyloseq.Rdata")


