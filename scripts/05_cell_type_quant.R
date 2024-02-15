library("tidyverse")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble") 
library("phyloseq")
library('vegan')
library('gridExtra')

path = "/home/user/Desktop/abhay/BIC_microbiome/results"
setwd(path)

custom_labeller <- function(variable, value) {
  labels <- c(
    "B_cells" = "B cells",
    "Cancer_fibroblasts" = "Cancer associated fibroblasts",
    "CD4_T_cells" = "CD4+ T-cells",
    "CD8_T_cells" = "CD8+ T-cells",
    "Endothlial_cells" = "Endothelial Cells",
    "Macrophages" = "Macrophages",
    "NK_cells" = "Natural Killer cells"
  )
  return(labels[value])
}

#Load EPIC immunescores
immune_scores <- read.csv("epic_immunedeconv_brca.csv", row.names = 1) 
immune_scores <- as.data.frame(t(immune_scores)) 
immune_scores$samples <- rownames(immune_scores) 
immune_scores <- inner_join(immune_scores, metadata, by = 'samples')


immune_scores <- immune_scores %>% 
  rename('CD4_T_cells' = `T cell CD4+`, 'CD8_T_cells' = `T cell CD8+`,
         'NK_cells' = `NK cell`, 'Macrophages' = `Macrophage`,
         'B_cells' = `B cell`, 'Cancer_fibroblasts' = `Cancer associated fibroblast`,
         'Endothlial_cells' = `Endothelial cell`) %>%
  select(samples, race, subtype, B_cells, Cancer_fibroblasts, CD4_T_cells, 
         CD8_T_cells, Endothlial_cells, Macrophages, NK_cells) %>%
  filter(subtype != 'NA') %>%
  filter(race != 'NA') %>% # count(race)
  pivot_longer(cols = c('B_cells', 'Cancer_fibroblasts', 'CD4_T_cells', 
                        'CD8_T_cells','Endothlial_cells', 'Macrophages', 'NK_cells'), 
               names_to = 'cell_type', values_to = 'cell_score') %>%
  mutate(race = factor(race, levels = c('Caucasian', 'Black', 'Asian'))) %>%
  mutate(subtype = factor(subtype, levels = c('Luminal-A', 'Luminal-B', 
                                              'Her2-pos', 'Triple-neg'))) %>%
  drop_na()

########################################################################
ggplot(immune_scores, aes(x=cell_type, y=cell_score, color = race)) +
  scale_color_manual(name=NULL,
                     values = c('blue', 'black', 'red'),
                     breaks = c('Caucasian', 'Black', 'Asian'),
                     labels = c('Caucasian (n = )', 'Black (n = )', 'Asian (n = )')) +
  scale_x_discrete(limits = c('B_cells', 'Cancer_fibroblasts', 'CD4_T_cells', 'CD8_T_cells', 'Endothlial_cells', 'Macrophages', 'NK_cells'),
                   labels = c('B cells', 'Cancer fibroblasts', 'CD4+ T-cells', 'CD8+ T-cells', 'Endothlial cells', 'Macrophages', 'NK cells')) +
  labs(title = NULL,
       x = 'Breast Cancer Immune Cell Type', 
       y = 'EPIC Cell Score') +
  theme_classic() +
  coord_cartesian(ylim=c(0,1))
########################################################################

# Significance testing 
immune_scores_wide <- immune_scores %>%
  pivot_wider(names_from = cell_type, values_from = cell_score) 

cell_types <- c('B_cells', 'Cancer_fibroblasts', 'CD4_T_cells', 'CD8_T_cells',
                'Endothlial_cells', 'Macrophages', 'NK_cells')

kt_results <- list()
pt_results <- list()

for (cell in cell_types) {
  kt_result <- kruskal.test(get(cell) ~ race, data = immune_scores_wide)
  
  if (kt_result$p.value < 0.05) {
    pt_result <- pairwise.wilcox.test(immune_scores_wide[[cell]],
                                      g = immune_scores_wide$race,
                                      p.adjust.method = 'BH')
    
    pt_results[[cell]] <- pt_result
  }
  
  kt_results[[cell]] <- kt_result
}

# Interpretation: Kruskal Walis (non-parametric), post-hoc: Pairwise Wilcox Test
# B cells             : Blacks and Caucasians (p = 1.424613e-05)
# Cancer fibroblasts  : Blacks and Caucasians (p = 3.137309e-04)
# CD4+ T cells        : ns
# CD8+ T cells        : Blacks and Caucasians (p = 0.018) 
#                       Blacks and Asians     (p = 0.018)
# Endothelial cells   : Blacks and Caucasians (p = 4.486494e-04)
# Macrophages         : Blacks and Caucasians (p = 4.486494e-04)
# NK cells            : Blacks and Caucasians (p = 7.627598e-04) 
#                       Caucasians and Asians (p = 8.9288232e-04)



# Facetted box-plot
 immune_scores %>%
  ggplot(aes(x = '', y = cell_score, fill = race)) +
  geom_boxplot() +
  scale_fill_manual(name=NULL,
                     values = c( 'gray','blue', 'red'),
                     breaks = c('Caucasian', 'Black', 'Asian'),
                     labels = c('Caucasian (n = 674)', 'Black (n = 158)', 'Asian (n = 59)')) +
  labs(
    title = NULL,
    x = 'Breast Cancer Immune Cell Type',
    y = 'EPIC Cell Score'
  ) +
  facet_wrap(~ cell_type, scales = "free", ncol = 4,
             labeller = labeller(cell_type = custom_labeller)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = c(0.87, 0.25),  
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 14),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.12)),
      scale_y_continuous(limits = c(0, 1.02)),
      scale_y_continuous(limits = c(0, 0.31)),
      scale_y_continuous(limits = c(0, 0.1)),
      scale_y_continuous(limits = c(0, 0.3)),
      scale_y_continuous(limits = c(0, 0.21)),
      scale_y_continuous(limits = c(0, 0.003)))
    )
ggsave('cell_quant.png', plot = last_plot(), dpi = 500, height = 8, width = 15)
  



