library("tidyverse")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble") 
library("phyloseq")
library('vegan')
library('gridExtra')

setwd("~/Desktop/AIIMS/BIC_microbiome")

######################################################
# Phyloseq installation and documentation
# library(BiocManager)
# BiocManager::install("phyloseq")
# browseVignettes("phyloseq") 
######################################################

#Load EPIC immunescores
immune_scores <- read.csv("cleaned_data/epic_immunedeconv_brca.csv", row.names = 1)
immune_scores<- as.data.frame(t(immune_scores)) 
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

# NMDS method for beta-diversity 
set.seed(1)
ndms_distance <- ordinate(physeq=bic, method="NMDS", distance = "bray")
ordinate.df <- as.data.frame(ndms_distance$points) %>%
  mutate(samples = rownames(.)) %>%
  inner_join(., metadata, by = 'samples') %>%
  filter(race != 'NA') %>%
  filter(race != "Native American") %>%
  filter(subtype != "NA")

ordinate.df %>%
  ggplot(aes(x = MDS1, y = MDS2, col = race, shape = subtype), size = 2) +
  geom_point(size = 2, alpha = 0.6) +
  labs(title = "Beta-diversity clustering by Race and Cancer subtype",
       x = 'NMDS1', y = 'NMDS2') +
  theme_classic()
ggsave("results/beta_diversity.png")


ordinate.df %>%
  ggplot(aes(x = MDS1, y = MDS2, col = race), size = 2) +
  geom_point() +
  labs(title = "Beta-diversity clustering by race",
       x = 'NMDS1', y = 'NMDS2') +
  facet_wrap(subtype~.) 
ggsave("results/beta_diversity_race.png")


ordinate.df %>%
  ggplot(aes(x = MDS1, y = MDS2, col = subtype), size = 2) +
  geom_point() +
  facet_wrap(race~.)
ggsave("results/beta_diversity_subtype.png")


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

p <- paste("P-value:", round(cor_shannon_CD4_T_cells$p.value, digits=5))
rho <- paste("rho:", round(cor_shannon_CD4_T_cells$estimate, digits=2))
annotation <- paste(p, rho, sep="\n")
ggplot(alpha_immune, aes(x=CD4_T_cells, y=Shannon)) +
  geom_point(aes(color=race)) +
  geom_smooth(method="lm", se=FALSE, color="black") +  # Use a single color for the line
  xlim(0, 0.20) +
  geom_text(aes(x=0.15, y=2, label=annotation), color="black", hjust = "left") +
  scale_color_manual(name=NULL,
                     values=c("grey","blue","red"),
                     breaks=c('Caucasian', "Black", "Asian"),
                     #labels=c('Caucasian', "Black", "Asian")
  ) +
  labs(title="",
       x="CD4+ T cell score (calculated by EPIC)",
       y="Shannon Diversity Index") +
  theme_classic()



plot_correlation <- function(cor_data, x_data, x_title, fig_title) {
  p <- paste("P-value:", round(cor_data$p.value, digits = 5))
  spearman_coeff <- paste("Spearman Coeff:", round(cor_data$estimate, digits = 2))
  annotation <- paste(p, spearman_coeff, sep = "\n")
  
  ggplot(alpha_immune, aes(x = x_data, y = Shannon)) +
    geom_point(aes(color = race)) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_text(aes(x = 0.15, y = 2, label = annotation),
              color = "black", hjust = "bottom", vjust = "top") +
    xlim(0, 0.20) +
    scale_color_manual(
      name = NULL,
      values = c("grey", "blue", "red"),
      breaks = c('Caucasian', "Black", "Asian")
    ) +
    labs(
      title = fig_title,
      x = x_title,
      y = "Shannon Diversity Index"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",  # Move the legend to the bottom
      legend.box.background = element_rect(colour = "black"),  # Add a black box around the legend
      legend.title = element_text(size = 14, face = "bold"),  # Increase legend title font size and make it bold
      plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins
      plot.title = element_text(hjust = 0.5)
    )
}



plot_correlation <- function(cor_data, x_data, x_title, fig_title) {
  p <- paste("P-value:", round(cor_data$p.value, digits = 5))
  spearman_coeff <- paste("Spearman Coeff:", round(cor_data$estimate, digits = 2))
  annotation <- paste(p, spearman_coeff, sep = "\n")
  
  ggplot(alpha_immune, aes(x = x_data, y = Shannon)) +
    geom_point(aes(color = race)) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_text(aes(x = 0.15, y = 2, label = annotation),
              color = "black", hjust = "bottom", vjust = "top", size = 6) +  # Increase annotation font size
    xlim(0, 0.25) +
    scale_color_manual(
      name = NULL,
      values = c("grey", "blue", "red"),
      breaks = c('Caucasian', "Black", "Asian")
    ) +
    labs(
      title = fig_title,  # Increase title font size and make it bold
      x = x_title,  # Use the provided x_title parameter for the x-axis label
      y = "Shannon Diversity Index"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",  # Move the legend to the bottom
      legend.box.background = element_rect(colour = "black"),  # Add a black box around the legend
      legend.title = element_text(size = 14, face = "bold"),  # Increase legend title font size and make it bold
      plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),  # Increase x-axis label font size
      axis.title.y = element_text(size = 14)  # Increase y-axis label font size
    )
}

CD4_race <- plot_correlation(cor_shannon_CD4_T_cells, 
                             alpha_immune$CD4_T_cells,
                             x_title = "Tumour infliitrating CD4+ T-cell score",
                             fig_title = "Correlation of Breast Cancer microbiome diversity\n with CD4+ T cell infiltration")

CD8_race <- plot_correlation(cor_shannon_CD8_T_cells, 
                             alpha_immune$CD8_T_cells,
                             x_title = "Tumour infliitrating CD8+ T-cell score",
                             fig_title = "Correlation of Breast Cancer microbiome diversity\n with CD8+ T cell infiltration")

# Reduce the size of the individual plots
CD4_race <- CD4_race + theme(plot.title = element_text(size = 16))  # Adjust title font size
CD8_race <- CD8_race + theme(plot.title = element_text(size = 16))  # Adjust title font size

T_cell_corr <- grid.arrange(CD4_race, CD8_race, ncol = 2)  # Reduce graph size
ggsave('shannon_Tcell_corr.png', T_cell_corr, height = 7, width = 15)


# write.csv(alpha_diversity, file = "results/brca_alpha_div.csv", col.names = T, row.names = T)
# save.image(file = "rdata/brca_phyloseq.Rdata")
load("rdata/brca_phyloseq.Rdata")


