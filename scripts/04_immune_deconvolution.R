library(immunedeconv)
library(pheatmap)
library(tidyverse)

path = "/home/user/Desktop/abhay/BIC_microbiome/"
setwd(path)

# Read RNAseq cnt_mtx for immune cell deconvolution
cnt_mtx <- read.csv("cleaned_data/brca_rnaseq.csv", row.names = 1)

# Method 1: TIMER
cancertype <- "BRCA"
res_timer = as.data.frame(deconvolute(cnt_mtx, "timer",indications=rep(tolower(cancertype),ncol(cnt_mtx))))
write.csv(res_timer, file = "results/timer_immunedeconv_brca.csv", row.names = F)

# Method 2: quanTIseq
res_quant = as.data.frame(deconvolute(cnt_mtx, "quantiseq"))
write.csv(res_quant, file = "results/quant_immunedeconv_brca.csv", row.names = F)

# Method 3: xCell
res_xcell = as.data.frame(deconvolute(cnt_mtx, "xcell"))
write.csv(res_xcell, file = "results/xcell_immunedeconv_brca.csv", row.names = F)

# Method 4: EPIC
res_epic = as.data.frame(deconvolute(cnt_mtx, "epic"))
write.csv(res_epic, file = "results/epic_immunedeconv_brca.csv", row.names = F)





