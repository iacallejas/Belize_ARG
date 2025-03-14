## CSV creation for Belize 2022 data by Ileana Galdamez (Callejas)
#Change this working drectory to the folder with the data and scripts
setwd("/Users/callejas/Library/CloudStorage/OneDrive-JPL/Belize 2022/Belize_ARG")
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

#IMPORT DATA
arg_data <- read_csv("Data/bz_arg_data.csv")

#NORMALIZE DATA BY 16SrRNA
arg_data <- arg_data %>%
  group_by(abv_name, full_name) %>%
  mutate(sul1_norm = sul1_abs/norm_16s,
         sul2_norm = sul2_abs/norm_16s,
         inti1_norm = inti1_abs/norm_16s,
         ermf_norm = ermf_abs/norm_16s,
         teta_norm = teta_abs/norm_16s,
         blashv_norm = blashv_abs/norm_16s) %>%
  ungroup()

#AVG PER SITE
arg_avg <-arg_data %>%
  group_by(abv_name, full_name) %>%
  summarize(sul1=mean(sul1_abs, na.rm = TRUE),
            sul2=mean(sul2_abs, na.rm = TRUE),
            inti1=mean(inti1_abs, na.rm = TRUE),
            ermf=mean(ermf_abs, na.rm = TRUE),
            teta=mean(teta_abs, na.rm = TRUE),
            blashv=mean(blashv_abs, na.rm = TRUE),
            norm_16s=mean(norm_16s, na.rm = TRUE),
            sul1_norm = mean(sul1_norm , na.rm = TRUE),
            sul2_norm = mean(sul2_norm , na.rm = TRUE),
            inti1_norm = mean(inti1_norm, na.rm = TRUE),
            ermf_norm = mean(ermf_norm, na.rm = TRUE),
            teta_norm = mean(teta_norm, na.rm = TRUE),
            blashv_norm = mean(blashv_norm, na.rm = TRUE)) %>%
  ungroup()

#Standard deviation #doesnt work for ones that only have one tube
arg_sd <- arg_data %>%
  group_by(abv_name, full_name) %>%
  summarize(sul1_sd=sd(sul1_abs, na.rm = TRUE),
            sul2_sd=sd(sul2_abs, na.rm = TRUE),
            inti1_sd=sd(inti1_abs, na.rm = TRUE),
            ermf_sd=sd(ermf_abs, na.rm = TRUE),
            teta_sd=sd(teta_abs, na.rm = TRUE),
            blashv_sd=sd(blashv_abs, na.rm = TRUE),
            norm_16s_sd=sd(norm_16s, na.rm = TRUE),
            sul1_norm_sd=sd(sul1_norm, na.rm = TRUE),
            sul2_norm_sd=sd(sul2_norm, na.rm = TRUE),
            inti1_norm_sd=sd(inti1_norm, na.rm = TRUE),
            ermf_norm_sd=sd(ermf_norm, na.rm = TRUE),
            teta_norm_sd=sd(teta_norm, na.rm = TRUE),
            blashv_norm_sd=sd(blashv_norm, na.rm = TRUE)) %>%
  ungroup()

#Export ARG AVG
#write.csv(arg_avg, "bz_arg_avg.csv")
