#Script for Metacompare plots by Ileana Galdamez (Callejas)
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
library(plot3D)
library(pals)

#IMPORT DATA
df <- read_csv("Data/MetaCompare.csv")

#ORDER DATA
df_ordered <- df %>%
  mutate(Site = factor(Site, 
                       levels = c("WCSR", "BRB", "FISH", "NA12", "NA11", "SNORKEL")))

# Create a "parula" color palette with 7 steps
col <- as.vector(ocean.phase(7))


# Resistance risk bar
ggplot(df_ordered, aes(x=Site, y=res_risk, fill=Site)) +
  geom_bar(stat="identity") +
  labs(title = "Resistome Risk Scores", y = "Resistome Risk\n", x = "Site")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank(),
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 35)) +
  scale_fill_manual(values=c("#A8780D","#CE2FCD","#D64759","#7E71F0","#339843","#1F93A7"))

## 3D Plot
scatter3D(x=df_ordered$`nARG/nContigs`, y=df_ordered$`nARG&MGE/nContigs` , z=df_ordered$`nARG&MGE&PAT/nContigs`, 
          phi = 5, theta = 40, bty = "g",
          pch = 16, cex = 2,
          col=c("#A8780D","#CE2FCD","#D64759","#7E71F0","#339843","#1F93A7"),
          colvar = as.integer(df_ordered$Site),
          main = "ARG, MGE, and Pathogen Data", xlab = "\nnARG/Contig",
          ylab ="\nnARG&MGE/Contig", zlab = "\nnARG&MGE&PATH/Contig")


text3D(x=df_ordered$`nARG/nContigs`, y=df_ordered$`nARG&MGE/nContigs` , z=df_ordered$`nARG&MGE&PAT/nContigs` ,
       labels = df_ordered$Site,
       add = TRUE, colkey = FALSE, cex = 1, col = "black")
