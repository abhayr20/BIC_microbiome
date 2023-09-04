# Load packages
library(immunedeconv)
library(maftools)
library(pheatmap)
library(tidyverse)

cnt_mtx <- read.csv("/data/abhay/brca/cleaned_data/brca_cnt_mtx.csv", row.names = 1)


quantiseq = as.data.frame(deconvolute(cnt_mtx, "quantiseq"))

# Method 1: CIBERSORT
res_ciber =

# Method 2: TIMER
cancertype <- "BRCA"
res_timer = as.data.frame(deconvolute(cnt_mtx, "timer",indications=rep(tolower(cancertype),ncol(cnt_mtx))))
write.csv(res_timer, file = "../results/timer_immunedeconv_brca.csv")

# Method 3: quanTIseq
res_quant = as.data.frame(deconvolute(cnt_mtx, "quantiseq"))
write.csv(res_quant, file = "../results/quant_immunedeconv_brca.csv")

# Method 4: xCell
res_xcell = as.data.frame(deconvolute(cnt_mtx, "xcell"))
write.csv(res_xcell, file = "../results/xcell_immunedeconv_brca.csv")

# Method 5: EPIC
res_epic = as.data.frame(deconvolute(cnt_mtx, "epic"))
write.csv(res_epic, file = "../results/epic_immunedeconv_brca.csv")

# Method 6: MPC-Counter
res_mcp = as.data.frame(deconvolute(cnt_mtx, "mcp_counter"))

