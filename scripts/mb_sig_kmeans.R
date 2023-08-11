library(tidyverse)
library(ggplot2)
library(factoextra)
library(gridExtra)

setwd("/Users/admin/Desktop/AIIMS/BIC_microbiome/")

# DATA PREPARATION
########################################################################
# Read raw data
mb <- read_tsv("data/BIC_microbiome_sig.tsv")
mb <- as.data.frame(mb)

# Rename the columns
mb_cols <- c("microbe", "mean_native_american", "mean_asian",
             "mean_black", "mean_white", "std_native_american", 
             "std_asian", "std_black", "std_white", 'p_value',
             "q_value", "enriched_population")
colnames(mb) <- mb_cols
rownames(mb) <- mb[,1]

# Filter significant data and remove irrelevant columns
mb_sig <- mb %>% 
  filter(p_value < 0.01) %>%
  select(mean_asian:mean_white) 

# compute variance of each variable
apply(mb_sig, 2, var)

#Scale the raw data
mb_sig <- scale(mb_sig)
head(mb_sig)

########################################################################

# PRINCIPAL COMPOONENT ANALYSIS
########################################################################

########################################################################

# K-MEANS CLUSTERING
########################################################################
#Evaluate clustering distance measures
distance <- get_dist(mb_sig)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

k2 <- kmeans(mb_sig, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2, data = mb_sig)

#Try to plot k = 3,4,5
k3 <- kmeans(mb_sig, centers = 3, nstart = 25)
k4 <- kmeans(mb_sig, centers = 4, nstart = 25)
k5 <- kmeans(mb_sig, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = mb_sig) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = mb_sig) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = mb_sig) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = mb_sig) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2)

# Determine optimal clusters 
fviz_nbclust(mb_sig, kmeans, method = "wss")          # Elbow plot
fviz_nbclust(mb_sig, kmeans, method = "silhouette")   # Silhoutte plot

########################################################################
