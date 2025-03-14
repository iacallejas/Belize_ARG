## Belize correlations by Ileana Galdamez (Callejas)
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
library(corrplot)
library(PerformanceAnalytics)
library(ggbiplot)
library(Hmisc)

#IMPORT DATA
df <- read_csv("Data/bz_all_data.csv")
type <- read_csv("Data/type_long.csv")
subtype <- read_csv("Data/subtype_long.csv")

#Metagenomics samples only
meta <- df %>%
  filter(abv_name %in% c("WCSR", "BRB", "FISH", "NA12", "NA11", "SNORKEL"))


#ONLY VALUES
#abs genes
numbers_abs <- df %>%
  select(toc, nitrate, phosphate, total_n, ammonia, ec, tc, esbl_ec, esbl_tc,
         ec_rr, sul1, sul2, inti1, ermf, teta, blashv, norm_16s)

colnames(numbers_abs) <- c("TOC", "Nitrate", "Phosphate", "Tot N", "Ammonia",
                           "E. coli", "Total Coliforms", "ESBL E. coli", "ESBL TC", "E. coli \n Resistance Ratio",
                           "sul1", "sul2", "intI1", "ermF", "tetA", "blaSHV", "16S rRNA")
cor_abs <- cor(numbers_abs, method = "pearson", use = "pairwise.complete.obs")

res2 <- rcorr(as.matrix(numbers_abs), type=c("pearson"))
res2$r
res2$P

res3=p.adjust(res2$P, method = "BH", 289)
res4=matrix(res3,nrow=17,ncol=17)

corrplot(res2$r, type = "upper", order = "hclust", p.mat = res2$P,
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F, 
         sig.level = c(.001, .01, .05), insig= "label_sig", tl.cex=1,
         pch.col = "yellow", pch.cex = 1)

#norm genes
numbers_norm <- df %>%
  select(toc, nitrate, phosphate, total_n, ammonia, ec, tc, esbl_ec, esbl_tc,
         ec_rr, norm_16s, sul1_norm, sul2_norm, inti1_norm, ermf_norm, teta_norm, blashv_norm)

colnames(numbers_norm) <- c("TOC", "Nitrate", "Phosphate", "Tot N", "Ammonia",
                           "E. coli", "Total Coliforms", "ESBL E. coli", "ESBL TC", "E. coli \n Resistance Ratio",
                           "16S rRNA", "sul1/16S", "sul2/16S", "intI1/16S", "ermF/16S", "tetA/16S", "blaSHV/16S")
cor_norm <- cor(numbers_norm, method = "pearson", use = "pairwise.complete.obs")


res2 <- rcorr(as.matrix(numbers_norm), type=c("pearson"))
res2$r
res2$P

res3=p.adjust(res2$P, method = "BH", 289)
res4=matrix(res3,nrow=17,ncol=17)

corrplot(res2$r, type = "upper", order = "hclust", p.mat = res2$P,
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F, 
         sig.level = c(.001, .01, .05), insig= "label_sig", tl.cex=1,
         pch.col = "yellow", pch.cex = 1)
#ABSOLUTE GENE COR PLOTS
#CORRPLOT W NUMBERS
corrplot(as.matrix(cor_abs), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black')

#WO diagonal
corrplot(as.matrix(cor_abs), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F)

#ANALYTICS
chart.Correlation(cor_abs, histogram=TRUE, pch=19)


#Normalized GENE COR PLOTS
#CORRPLOT W NUMBERS
corrplot(as.matrix(cor_norm), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black')

#WO diagonal
corrplot(as.matrix(cor_norm), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F)

#ANALYTICS
chart.Correlation(cor_norm, histogram=TRUE, pch=19)


##############
## Removing FIB due to low values

numbers_abs <- df %>%
  select(toc, nitrate, phosphate, total_n, ammonia,
         sul1, sul2, inti1, ermf, teta, blashv, norm_16s)

colnames(numbers_abs) <- c("TOC", "Nitrate", "Phosphate", "Tot N", "Ammonia",
                           "sul1", "sul2", "intI1", "ermF", "tetA", "blaSHV", "16S rRNA")

cor_abs <- cor(numbers_abs, method = "pearson", use = "pairwise.complete.obs")

#norm genes
numbers_norm <- df %>%
  select(toc, nitrate, phosphate, total_n, ammonia,
         norm_16s, sul1_norm, sul2_norm, inti1_norm, ermf_norm, teta_norm, blashv_norm)

colnames(numbers_norm) <- c("TOC", "Nitrate", "Phosphate", "Tot N", "Ammonia",
                            "16S rRNA", "sul1/16S", "sul2/16S", "intI1/16S", "ermF/16S", "tetA/16S", "blaSHV/16S")
cor_norm <- cor(numbers_norm, method = "pearson", use = "pairwise.complete.obs")

#ABSOLUTE GENE COR PLOTS
#CORRPLOT W NUMBERS
corrplot(as.matrix(cor_abs), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black')

#WO diagonal
corrplot(as.matrix(cor_abs), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F)

#ANALYTICS
chart.Correlation(cor_abs, histogram=TRUE, pch=19)


#Normalized GENE COR PLOTS
#CORRPLOT W NUMBERS
corrplot(as.matrix(cor_norm), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black')

#WO diagonal
corrplot(as.matrix(cor_norm), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,addCoef.col = 'black', diag=F)

#ANALYTICS
chart.Correlation(cor_norm, histogram=TRUE, pch=19)


## PCA analysis
#All genes
pca <- prcomp(df[,c(4,5,8,18:30)], center=TRUE, scale=TRUE)

summary(pca)

ggbiplot(pca, ellipse=F, labels=df$abv_name, groups = df$group) +
  theme_light()


#Abs genes
pca <- prcomp(df[,c(4,5,8,18:24)], center=TRUE, scale=TRUE)

summary(pca)

ggbiplot(pca, ellipse=F, labels=df$abv_name, groups = df$group) +
  theme_light()

#Norm genes
pca <- prcomp(df[,c(4,5,8,25:30)], center=TRUE, scale=TRUE)

summary(pca)

ggbiplot(pca, ellipse=F, labels=df$abv_name, groups = df$group) +
  theme_light()

