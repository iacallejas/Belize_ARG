#Belize Heatplots between methods by Ileana Galdamez (Callejas)
setwd("/Users/callejas/Library/CloudStorage/OneDrive-JPL/Belize\ 2022/Belize_ARG")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(plotly)
library(scales)
library(ggpubr)
library(plotrix)
library(stringr)
library(pheatmap)
library(pals)


#IMPORT DATA
arg_avg <- read_csv("Data/bz_arg_avg.csv")

#SELECT NECESSARY COLUMNS
qpcr <- arg_avg %>%
  select("abv_name","full_name", "sul1_norm", "sul2_norm",
         "inti1_norm","ermf_norm","teta_norm","blashv_norm")

#FILTER BY NAME
qpcr_sub <- qpcr %>%
  filter(abv_name %in% c("NA11", "NA12", "POUL","FISH", "HCM",
                         "MAN", "BRB", "BRO","WCSR","WCSL",
                         "SNORKEL"))%>%
  dplyr::rename(sul1 = sul1_norm,
                sul2 = sul2_norm,
                intI1 = inti1_norm,
                ermF = ermf_norm,
                tetA = teta_norm,
                blaSHV = blashv_norm)
#CREATE MATRIX
qpcr_mat <- qpcr_sub %>%
  select(3:8)
rownames(qpcr_mat) <- paste(qpcr_sub$abv_name)

#REORDER
qpcr_mat <- qpcr_mat %>%
  slice(6,7,8,3,4,5,1,2,11,10,9)

#METAGENOMICS HEATMAP
col <- as.vector(rev(ocean.deep(100)))
pheatmap(log10(as.matrix(qpcr_mat)), labels_row= rownames(qpcr_mat), 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         color = col, angle_col = c("0"),
         main = "qPCR")

####################################################################
########### METGENOMICS
####################################################################
#IMPORT METAGENOMICS DATA
meta_long <- read_csv("Data/type_long.csv")

meta_sub <- meta_long %>%
  pivot_wider(names_from = `AMR Class`, values_from = norm16s)

#CREATE MATRIX
meta_mat <- meta_sub %>%
  select(2:18)
rownames(meta_mat) <- paste(meta_sub$name)

#REORDER
meta_mat <- meta_mat %>%
  slice(6,3,1,4,2,5)

# Replace 1e-10 with an appropriate small constant based on the scale of your data
small_const <- 0.0000001

# Adding the small constant to all elements of fibdata
meta_mat_adj <- meta_mat + small_const

#QPCR HEATMAP
col2 <- as.vector(rev(magma(100)))
pheatmap(log10(as.matrix(meta_mat_adj)), labels_row= rownames(meta_mat_adj), 
         cluster_rows=FALSE, cluster_cols=FALSE, color = col2,
         angle_col = c("90"),
         main = "Metagenomics")

##################################################################
########## Culture
##################################################################
#IMPORT DATA
idexx <- read_csv("Data/idexx.csv")

#SELECT NECESSARY COLUMNS
cul_sub <- idexx %>%
  select("Site","ESBL-EC")

#CREATE MATRIX
cul_mat <- cul_sub %>%
  select(2)
rownames(cul_mat) <- paste(cul_sub$Site)

#REORDER
cul_mat <- cul_mat %>%
  slice(2,1,5,6,4,3,7,8,10,9)

# Replace 1e-10 with an appropriate small constant based on the scale of your data
small_const2 <- 1

# Adding the small constant to all elements of fibdata
cul_mat_adj <- cul_mat + small_const2

#ESBL ECOLI HEATMAP
col3 <- as.vector(rev(ocean.ice(100)))
pheatmap(log10(as.matrix(cul_mat_adj)), labels_row= rownames(cul_mat_adj), 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         color = col3, angle_col = c("0"),
         main = "Culture")

