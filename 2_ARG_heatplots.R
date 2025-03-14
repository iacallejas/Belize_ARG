## Belize ARG Trends by Ileana Galdamez (Callejas)
#Set working directory
setwd("/Users/callejas/Library/CloudStorage/OneDrive-JPL/Belize 2022/Belize_ARG ")
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
library(grid)

#IMPORT DATA
arg_avg <- read_csv("Data/bz_arg_avg.csv")

arg_avg <- arg_avg %>%
  slice(16,18,1,2,8,6,12,4,10,9,11,5,15,17,3,7,13,14)

bz_arg_data <- read_csv("bz_arg_data.csv")

#CREATE MATRIX
#ABS
data_abs <- arg_avg  %>%
  select(c('full_name', 'sul1', 'sul2', 'inti1', 'ermf', 'teta', 'blashv')) %>%
  select(2:7)


rownames(data_abs) <- paste(arg_avg$full_name)

#NORM
data_norm <- arg_avg  %>%
  select(c('full_name', 'sul1_norm', 'sul2_norm', 'inti1_norm', 'ermf_norm', 'teta_norm', 'blashv_norm')) %>%
  select(2:7)

rownames(data_norm) <- paste(arg_avg$full_name)

#HEATPLOTS
#ABS **Updated
col <- as.vector(ocean.balance(100))
colnames(data_abs) <- c("sul1", "sul2", "intI1", "ermF", "tetA", "blaSHV")
pheatmap(log10(data_abs), labels_row= rownames(data_abs), angle_col = 0, color = col)
grid.text("log(GC/mL)", x = 0.95, y = 0.9, rot = 0, gp = gpar(fontsize = 10))

#NORM
pheatmap(log10(data_norm), labels_row= rownames(data_norm))
