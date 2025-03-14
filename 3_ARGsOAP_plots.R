## BZ ARGsOAP data by Ileana Galdamez (Callejas)
setwd("/Users/callejas/Library/CloudStorage/OneDrive-JPL/Belize 2022/Belize_ARG")
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(pheatmap)
library(stringr)
library(pals)
library(circlize)

# READ IN FILES
subtype <- read_tsv("Data/extracted-1.fa.normalize_16s.subtype.tab.txt")
type <- read_tsv("Data/extracted-1.fa.normalize_16s.type.tab.txt")

#RENAME COLUMNS
subtype <- subtype %>%
  dplyr::rename("Gene" = "Subtype_to_16S_copy_number",
      "FISH" = "A1",
         "WCSR" = "A2",
         "NA12" = "A3",
         "BRB" = "A4",
         "SNORKEL" = "A5",
         "NA11" = "A6")

type <- type %>%
  dplyr::rename("AMR Class" = "Type_to_16S_copy_number",
    "FISH" = "A1",
                "WCSR" = "A2",
                "NA12" = "A3",
                "BRB" = "A4",
                "SNORKEL" = "A5",
                "NA11" = "A6")

# CAPITALIZE ALL ANTIBIOTICS
type$`AMR Class` <- str_to_title(type$`AMR Class`)

#RENAME MLS
type <- type %>%
  dplyr::mutate(`AMR Class` = ifelse(`AMR Class` == "Macrolide-Lincosamide-Streptogramin", "MLS", `AMR Class`))


#REMOVE ZERO DATA
subtype_clean <- subtype[rowSums(subtype[,2:7])>0,]

type_clean <- type[rowSums(type[,2:7])>0,]


#STACKED BAR PLOTS
#Change data to long format
type_long <- gather(type_clean, name, norm16s, -`AMR Class`)
subtype_long <- gather(subtype_clean, name, norm16s, -`Gene`)
#write.csv(type_long, "type_long.csv", row.names = FALSE)
#write.csv(subtype_long, "subtype_long.csv", row.names = FALSE)


#Reorder names to match metacompare
type_long <- type_long %>%
  mutate(name = factor(name, 
                       levels = c("WCSR", "BRB", "FISH",
                                  "NA12", "NA11", "SNORKEL")))


#Stacked bar plot
ggplot(type_long, aes(fill=`AMR Class`, y=norm16s, x=name)) + 
  geom_bar(stat="identity") +
  labs(y = "Relative Abundance", x = "Site", fill="AMR Class")+
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key=element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, .2)) +
  scale_fill_manual(values=as.vector(stepped3(17)))


###########
# CHORD PLOT

#Set groups
type_long <- type_long %>%
  mutate(Group = case_when(
    name %in% c("WCSR", "BRB") ~ "Belize River",
    name %in% c("FISH", "NA12") ~ "Haulover Creek in Belize City",
    name == "SNORKEL" ~ "Coral Reef",
    name == "NA11" ~ "Treatment Lagoon",
    TRUE ~ NA_character_
  ))

# CALC AVERAGES PER LAND TYPE
type_avg_group <- type_long %>%
  dplyr::group_by(`AMR Class`, Group)%>%
  dplyr::summarize(mean_16s = mean(norm16s, na.rm = TRUE))%>%
  ungroup()

# Find top 10 AMR Class by mean_16s
top10_AMR_class <- type_avg_group %>%
  dplyr::group_by(`AMR Class`) %>%
  dplyr::summarize(sum=sum(mean_16s)) %>% 
  top_n(10,sum) %>%
  arrange(desc(sum))

print(top10_AMR_class)

#Filter based on top10
filtered_df <- type_long %>%
  filter(`AMR Class` %in% top10_AMR_class$`AMR Class`) %>%
  select(-name) %>%
  select(`AMR Class`,Group,norm16s)

#######
# TOP 10 SUBTYPE
# CALC AVERAGES PER LAND TYPE
subtype_avg_group <- subtype_long %>%
  dplyr::group_by(`Gene`, name)%>%
  dplyr::summarize(mean_16s = mean(norm16s, na.rm = TRUE))%>%
  ungroup()

# Find top 10 AMR Class by mean_16s
top10_AMR_class_subtype <- subtype_avg_group %>%
  dplyr::group_by(`Gene`) %>%
  dplyr::summarize(sum=sum(mean_16s)) %>% 
  top_n(10,sum) %>%
  arrange(desc(sum))

####################################
# AVERAGE PER SITE
type_avg_name <- type_long %>%
  dplyr::group_by(`AMR Class`, name)%>%
  dplyr::summarize(mean_16s = mean(norm16s, na.rm = TRUE))%>%
  ungroup()

# Find top 10 AMR Class by mean_16s
top10_AMR_class_name <- type_avg_name %>%
  dplyr::group_by(`AMR Class`) %>%
  dplyr::summarize(sum=sum(mean_16s)) %>% 
  top_n(10,sum) %>%
  arrange(desc(sum))

print(top10_AMR_class_name)

#Filter based on top10
filtered_df_name <- type_long %>%
  filter(`AMR Class` %in% top10_AMR_class_name$`AMR Class`) %>%
  select(-Group) %>%
  select(`AMR Class`,name,norm16s)

#Chord plot
colors <- c(WCSR = "#A8780D", BRB = "#CE2FCD",
            FISH = "#D64759", NA12= "#7E71F0",
            NA11 = "#339843", SNORKEL="#1F93A7",
            Bacitracin = "#440154",
            Multidrug = "#482878", Unclassified = "#3e4989", 
            MLS = "#31688e", Tetracycline = "#26828e",
            `Beta-Lactam` = "#1f9e89",  Sulfonamide= "#35b779", 
            Vancomycin = "#6ece58", Aminoglycoside = "#b5de2b", Rifamycin = "#fde725")
circos.par(start.degree = 90)
chordDiagram(filtered_df_name, 
             order = c("WCSR", "BRB",
                       "FISH", "NA12", "NA11", "SNORKEL",
                       "Bacitracin","Multidrug","Unclassified","MLS",           
                       "Tetracycline","Beta-Lactam","Sulfonamide","Vancomycin",
                       "Aminoglycoside","Rifamycin"),
             grid.col = colors, transparency = 0.25,
             big.gap = 20, annotationTrack = "grid", preAllocateTracks = 1, 
             directional = 1,
             direction.type = c("arrows", "diffHeight"), 
             diffHeight  = -0.01,
             link.arr.type = "big.arrow",
             link.sort = TRUE, link.decreasing = FALSE, link.largest.ontop = TRUE)


circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .3, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.6, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()





