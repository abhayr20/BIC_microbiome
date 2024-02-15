library(tidyverse)
library(gtsummary)
library(ggtext)
library(RColorBrewer)
library(vegan)

setwd("path/to/raw_data")

brca_subtype <- data.frame(
  subtype = c('Luminal-A', 'Luminal-B', 'Her2-pos', 'Triple-neg'),
  ER_and_PR =  c("Atleast 1 pos","Atleast 1 pos", "neg", "neg"),
  HER = c("neg", "pos/neg", "pos", "neg"),
  Prognosis = c('good', "mid", "mid/bad", "poor")
) 

brca_map <- read_csv("raw_data/brca_map.csv") %>%
  select(X, subtype) %>%
  rename(samples = X)

metadata <- read_csv("raw_data/metadata.csv") %>%
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

meta_race <- metadata %>%
  filter(!race %in% c('NA', "Native American")) 

meta_subtype <- metadata %>%
  filter(!subtype %in% c('NA')) 

taxonomy <- read_csv("raw_data/taxonomy.csv")

otu_counts <- read_csv("raw_data/otu.csv") %>%
  rename(otu = 1) %>%
  pivot_longer(cols = -otu, names_to = 'samples') %>% 
  pivot_wider(names_from = otu, values_from = value) %>%
  pivot_longer(-samples) %>%
  group_by(samples) %>%
  mutate(total = sum(value)) %>%
  filter(total > 100000) %>%       # Removes 10 samples with less than 100k OTUs
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%           # Removes 165 OTUs with zero reads
  ungroup() %>%
  select(-total) %>%
  pivot_wider(samples) %>%
  as.data.frame()
  
rownames(otu_counts) <- otu_counts$samples
otu_counts <- otu_counts[, -1]
otu_counts <- as.matrix(otu_counts) 

#Rarefying OTU counts 
#Controls uneven sampling effort by sub sampling 100000 OTUs 100 times
set.seed(19970911)
dist <-  avgdist(otu_counts, demethod = "bray", sample = 100000) 
#dist <- vegdist(otu_counts, method = "bray")

# Calculate ordination
set.seed(4)
nmds <- metaMDS(dist)

scores(nmds) %>%
  as_tibble(rownames = "samples") %>%
  inner_join(., metadata, by = 'samples') %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=subtype)) +
  geom_point()

save.image(file = 'rdata/rarefied.Rdata')
#load('rdata/rarefied.Rdata')

# Alpha Diversity
otu_counts

