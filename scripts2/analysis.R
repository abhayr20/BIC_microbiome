library(maftools)
library(pheatmap)
library(tidyverse)
library(corrplot)

setwd("/data/abhay/brca")

#read data
xcell     <- read.csv("results/xcell_immunedeconv_brca.csv")
quantiseq <- read.csv("results/quant_immunedeconv_brca.csv")
timer     <- read.csv("results/timer_immunedeconv_brca.csv")
epic      <- read.csv("results/epic_immunedeconv_brca.csv")
alpha     <- read.csv("results/brca_alpha_div.csv", row.names = 1)
brca_map  <- read.csv("cleaned_data/brca_map.csv")

perform_analysis <- function(df, alpha) {
  # Extract the name of the dataframe
  df_name <- deparse(substitute(df))
  
  # Transpose the input dataframe
  df <- t(df)
  colnames(df) <- unlist(df[1, ])
  df <- df[-1, ]
  df <- as.data.frame(df)
  
  # Match row names between df and alpha
  alpha_sort <- alpha[match(rownames(df), rownames(alpha)), ]
  
  # Check if row names match
  if (!identical(rownames(alpha_sort), rownames(df))) {
    stop("Row names do not match between df and alpha.")
  }
  
  # Add the alpha_div column from alpha to df
  df$alpha_div <- alpha_sort$Shannon
  
  # Convert columns to numeric
  for (col in names(df)) {
    df[[col]] <- as.numeric(df[[col]])
  }
  
  # Calculate the Spearman correlation matrix
  corr_df <- cor(df, method = 'spearman')
  
  # Create a correlation plot with a dynamic title
  title <- paste(df_name, ": Correlation Plot Of Immune cell signatures with alpha-diversity", sep="")
  corrplot(corr_df, tl.col = "brown", tl.srt = 30, bg = "White",
           title = title, type = "full")
  
  # Save the correlation plot as a PNG file
  png(file = paste0("results/", df_name, "_immunecell_corr.png"))
  corrplot(corr_df, tl.col = "brown", tl.srt = 30, bg = "White",
           title = title, type = "full")
  dev.off()
  
  return(df)
}

perform_analysis(epic, alpha)
perform_analysis(xcell, alpha)
perform_analysis(timer, alpha)
perform_analysis(quantiseq, alpha)

ggplot(brca_map, aes(x = Race.category)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", na.rm = TRUE) +
  labs(title = "Histogram of Values (Ignoring NAs)", x = "Values", y = "Frequency")

1. Rename function to global correlation
2. Make the correlation plot