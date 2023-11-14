library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/Biotin Deficiency 2022 Final/")

### Trios ---
metadata <- read.table("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
pathway <- read.delim("Trios/differential_Pathway/feature-table.tsv", header = TRUE, row.names=1)
pathway <- pathway %>% select(-c("taxonomy"))

# Get datasets ---
input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter((Subset =="Mucosal_Colon" & Genotype !="HET"), SampleID %in% names(df_input_data)) %>% pull(SampleID)
df_input_data <- pathway[, samples]

heme_data <- as.data.frame(t(df_input_data))
heme_data <- heme_data %>% select(c("PWY-5920"))
heme_data$SampleID <- row.names(heme_data)
heme_data <- merge(heme_data,input_metadata, by="SampleID")

cols <- c("MUT"="red", "WT"="black")

ggplot(heme_data, aes(`PWY-5920`, fill = Genotype, colour = Genotype)) +
  geom_density(alpha = 0.5) + 
  theme_cowplot(16) +
  scale_x_log10()+ 
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  labs(title = "Heme biosynthesis",
       x = "log10(Abundance)") + 
  theme_cowplot(16) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "top", legend.justification = "center")+
  theme(legend.title=element_blank())

# Make a density plot for STool ---
cecum_data <- as.data.frame(t(stool_data))
cecum_data <- cecum_data %>% select(c("BIOTIN-BIOSYNTHESIS-PWY"))
cecum_data$SampleID <- row.names(cecum_data)
cecum_data <- merge(cecum_data,input_metadata, by="SampleID")
cecum_data$Group<-revalue(cecum_data$Group, c("Control"="CD", "BD"="BD"))
cols <- c("BD"="red", "CD"="black")

ucla_stool_biotin <- ggplot(cecum_data, aes(`BIOTIN-BIOSYNTHESIS-PWY`, fill = Group, colour = Group)) +
  geom_density(alpha = 0.5) + 
  theme_cowplot(16) +
  scale_x_log10()+ 
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  labs(title = "Biotin biosynthesis I",
       x = "log10(Abundance)") + 
  theme_cowplot(16) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "top", legend.justification = "center")+
  theme(legend.title=element_blank())


### UCI grab biotin biosynthesis ---
input_data <- read.csv("UCI/pathway-table/UCI Metacyc Table - uci_pwy.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
annotation <- input_data %>% select("description")
annotation$feature <- row.names(input_data)
df_input_data <- select(input_data, -c("description", "taxonomy"))

input_metadata <-read.delim("UCI/starting_files/UCI_metadata_analysis_nofood.tsv",sep="\t",header=TRUE, row.names=1)
input_metadata$SampleID <- row.names(input_metadata)

# Get datasets ---
samples <- input_metadata %>% filter(Sample_type =="stool", SampleID %in% names(df_input_data)) %>% pull(SampleID)
stool_data <- df_input_data[,samples]

samples <- input_metadata %>% filter(Sample_type =="intestine", SampleID %in% names(df_input_data)) %>% pull(SampleID)
colon_data <- df_input_data[,samples]

# Make a density plot ---
colon_data <- as.data.frame(t(colon_data))
colon_data <- colon_data %>% select(c("BIOTIN-BIOSYNTHESIS-PWY"))
colon_data$SampleID <- row.names(colon_data)
colon_data <- merge(colon_data,input_metadata, by="SampleID")

cols <- c("KO"="red", "WT"="black")

uci_colon_biotin <- ggplot(colon_data, aes(`BIOTIN-BIOSYNTHESIS-PWY`, fill = Genotype, colour = Genotype)) +
  geom_density(alpha = 0.5) + 
  theme_cowplot(16) +
  scale_x_log10()+ 
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  labs(title = "Biotin biosynthesis I",
       x = "log10(Abundance)") + 
  theme_cowplot(16) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "top", legend.justification = "center")+
  theme(legend.title=element_blank())

first <- plot_grid(uci_colon_pwy, ucla_cecum_pwy, ucla_stool_pwy, ucla_stool_biotin, nrow=3, labels=c("A","B", "C"))
second <- plot_grid(uci_colon_biotin, ucla_cecum_biotin, NULL, nrow=3,align="hv", axis="tblr", labels=c("D", "E", ""))
dev.new(width=10, height=20)
plot_grid(first, second, axis="tblr", align="hv")

third<- plot_grid(ucla_cecum_pwy, ucla_cecum_biotin, 
                  ucla_stool_pwy, ucla_stool_biotin, 
                  uci_colon_pwy, uci_colon_biotin, nrow=3, labels=c("A","D","B","E","C","F"))
fourth <- plot_grid(uci_colon_pwy, uci_colon_biotin,  nrow=1, labels=c("C",""))

biotin_only <- plot_grid(ucla_cecum_biotin, 
                        ucla_stool_biotin, 
                        uci_colon_biotin, ncol=3)

dev.new(width=10, height=20)
third

dev.new(width=10, height=20)
fourth