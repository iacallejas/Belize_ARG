#FIB plots by Ileana Galdamez (Callejas)
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
library(cowplot)

#IMPORT DATA
fib <- read_csv("Data/idexx.csv")

fib <- fib %>%
  mutate(per_esbl_rr = `EC-RR`*100)

fib_ordered <- fib %>%
  mutate(Site = factor(Site, 
                       levels = c("BRB", "BRO", "MAN", "HCM", "POUL","FISH","NA12","NA11")))


################################
################################
################# Plots
################################
################################

## TC
tc_plot <- ggplot(fib_ordered, aes(x = Site, y=TC)) +
  geom_bar(stat='identity', fill = "dodgerblue") +
  geom_hline(yintercept = 10000, col = "red") +
  labs(title = "Total Coliform Levels",  y = "log MPN/100 mL", x = "Site")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits=c(1, 100000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0,0))

tc_plot

## EC
ec_plot <- ggplot(fib_ordered, aes(x = Site, y=EC)) +
  geom_bar(stat='identity', fill = "dodgerblue") +
  geom_hline(yintercept = 235, col = "yellow") +
  geom_hline(yintercept = 576, col = "red") +
  labs(title = "E. coli Levels",  y = "log MPN/100 mL", x = "Site")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits=c(1, 100000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0,0))

ec_plot

## ESBL EC
esbl_ec_plot <- ggplot(fib_ordered, aes(x = Site, y=`ESBL-EC`)) +
  geom_bar(stat='identity', fill = "dodgerblue") +
  labs(title = "ESBL E. coli Levels",  y = "MPN/100 mL", x = "Site")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) 

esbl_ec_plot

## ESBL EC RR
ec_rr_plot <- ggplot(fib_ordered, aes(x = Site, y=per_esbl_rr)) +
  geom_bar(stat="identity", fill = "dodgerblue") +
  labs(title = "ESBL E. coli Resistance Ratio",  y = "Resistance Ratio (%)", x = "Site")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank(),
        legend.title = element_blank())  + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 10)) 

ec_rr_plot

# arrange plots in a 2x2 grid
plot_grid(tc_plot, ec_plot, esbl_ec_plot, ec_rr_plot, ncol = 2)

