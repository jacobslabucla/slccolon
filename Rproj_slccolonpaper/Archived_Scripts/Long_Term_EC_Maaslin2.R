library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Long_Term_PWY_Maaslin2.R")


metadata <- read.csv("Long_Term/starting_files/SLC_LT_metadata.csv", header=TRUE)
counts <- read.table("Long_Term/starting_files/SLT_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)
pathway <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_ec_metagenome/feature-table.tsv", header = TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Grab samples in the pathway table that are present in counts --
names(pathway)
pathway <- pathway %>% select(c(names(counts)))

### Luminal Colon ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Subset =="Luminal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon", "Cecum"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID"), 
                    fixed_effects = c("Site",  "Sex", "Genotype"), normalization = "TSS", 
                    random_effects = c("MouseID"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)


### Mucosal Colon ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Subset =="Mucosal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon", "Cecum"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Long_Term/differential_EC/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID"), 
                    fixed_effects = c("Site",  "Sex", "Genotype"), normalization = "TSS", 
                    random_effects = c("MouseID"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)



### Visualize PWY results ---

## Luminal Colon --
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_ec_metagenome/annotated_ec.tsv",row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
data <- merge(data,annotation, by="feature")
write.csv(data, "Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/annotated_significant_results_pwy_q0.25.tsv")

res_plot <- data %>% filter(value=="MUT")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "MUT"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
pwy_mut <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.25, abs(coef) > 0.5) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (MUT/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("Luminal Colon: WT vs MUT q<0.25") +
  theme(plot.title = element_text(hjust = 0.5))
pwy_mut

res_plot <- data %>% filter(value=="HET")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "HET"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
pwy_het <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.25, abs(coef) > 0.5) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (HET/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("Luminal Colon: WT vs HET q<0.25") +
  theme(plot.title = element_text(hjust = 0.5))
pwy_het

plot_grid(pwy_mut, pwy_het, labels=c("A","B"))


## Mucosal Colon --
data<-read.table("Long_Term/differential_EC/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_ec_metagenome/annotated_ec.tsv",row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
data <- merge(data,annotation, by="feature")
write.csv(data, "Long_Term/differential_EC/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/annotated_significant_results_pwy_q0.25.tsv")

res_plot <- data %>% filter(value=="MUT")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "MUT"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
mc_pwy_mut <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.25, abs(coef) > 0.5) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (MUT/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("Mucosal Colon: WT vs MUT q<0.25") +
  theme(plot.title = element_text(hjust = 0.5))

res_plot <- data %>% filter(value=="HET")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "HET"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
mc_pwy_het <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.25, abs(coef) > 0.5) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (HET/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("Mucosal Colon: WT vs HET q<0.25") +
  theme(plot.title = element_text(hjust = 0.5))
pwy_het

plot_grid(mc_pwy_mut, mc_pwy_het, labels=c("A","B"))
