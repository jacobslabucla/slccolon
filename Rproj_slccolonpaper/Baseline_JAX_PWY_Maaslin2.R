library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Baseline_JAX_PWY_Maaslin2.R")

metadata <- read.table("Baseline/starting_files/Baseline_Metadata.tsv", header=TRUE)
counts <- read.table("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)
pathway <- read.delim("Baseline/differential_Pathway/feature-table.tsv", header = TRUE, row.names=1)
enzyme <- read.delim("Baseline/differential_Pathway/feature-table.tsv", header = TRUE, row.names=1)
ko <- read.delim("Baseline/differential_KO/feature-table.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Get counts only from the JAX mice

JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts))
row.names(JAX_meta) <- JAX_meta$SampleID
JAX <- JAX_meta$SampleID
JAX_counts <- counts %>% select(all_of(JAX))

## Grab samples in the pathway table that are present in counts
names(pathway)
pathway <- pathway %>% select(c(names(JAX_counts)))
enzyme <- enzyme %>% select(c(names(JAX_counts)))
ko <-ko %>% select(c(names(JAX_counts)))

### Run Maaslin2 and get table of relative abundances

run_Maaslin2 <- function(counts, metadata, subset_string) {
  #input_data <- read.delim("export_s3_min10000_PFF_Baseline_min10000_no_tax_PFF_ASV_table/feature-table.tsv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
  input_data <- as.data.frame(counts)
  df_input_data <- as.data.frame(input_data)

  input_metadata <- as.data.frame(metadata)

  target <- colnames(df_input_data)
  input_metadata = input_metadata[match(target, row.names(input_metadata)),]
  target == row.names(input_metadata)
  
  df_input_metadata<-input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
  df_input_metadata$Sex <- factor(df_input_metadata$Sex)
  sapply(df_input_metadata,levels)
  
  ?Maaslin2
  fit_data = Maaslin2(input_data=df_input_data, 
                      input_metadata=df_input_metadata, 
                      output = paste0("Baseline/",subset_string,"_Maaslin2_Sex_Genotype"), 
                      fixed_effects = c("Sex","Genotype"),normalization="TSS", 
                      min_prevalence = 0.15,
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}


## Sex and Genotype

# Baseline
run_Maaslin2(counts = pathway,
             metadata= JAX_meta,
             subset_string = "Baseline_JAX_PWY")
run_Maaslin2(counts = enzyme,
             metadata= JAX_meta,
             subset_string = "Baseline_JAX_EC")
run_Maaslin2(counts = ko,
             metadata= JAX_meta,
             subset_string = "Baseline_JAX_KO")

### Visualize EC results  ---

# Baseline WT vs HET
data<-read.table("Baseline/differential_Pathway/feature-table.tsv", header=TRUE)
data <- data %>% filter(qval <0.10)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_ec_metagenome/annotated_EC.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>% filter(value=="HET")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "HET"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
baseline_ec_het <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.10, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (HET/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("WT vs HET: EC data") +
  theme(plot.title = element_text(hjust = 0.5))
baseline_ec_het

# Baseline WT vs MUT
data<-read.table("ASV-level_SLC_Baseline_EC_Maaslin2_Line_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.10)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_ec_metagenome/annotated_EC.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>% filter(value=="MUT")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "MUT"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
baseline_ec_mut <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.10, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (MUT/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("WT vs MUT: EC data") +
  theme(plot.title = element_text(hjust = 0.5))
baseline_ec_mut

plot_grid(baseline_ec_het, baseline_ec_mut, labels = c("A", "B"))
### Visualize Pathway results  ---

# Baseline WT vs HET
data<-read.table("ASV-level_SLC_Baseline_PWY_Maaslin2_Line_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.20)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_pathway_abundance/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>% filter(value=="HET") 
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "HET"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
baseline_pwy_het <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.20, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (HET/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("WT vs HET: MetaCyc data") +
  theme(plot.title = element_text(hjust = 0.5))
baseline_pwy_het


# Baseline WT vs MUT
data<-read.table("ASV-level_SLC_Baseline_PWY_Maaslin2_Line_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.20)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_pathway_abundance/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)

data <- merge(data,annotation, by="feature")

res_plot <- data %>% filter(value=="MUT")
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "WT", "MUT"))

y = tapply(res_plot$coef, res_plot$description, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$description= factor(as.character(res_plot$description), levels = names(y))

cols <- c("WT"="black", "HET"="blue", "MUT"="firebrick")
baseline_pwy_mut <- res_plot %>%
  arrange(coef) %>%
  filter(qval < 0.20, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, description, fill = site)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(16) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (MUT/WT)",
       y = "",
       fill = "") +
  theme(legend.position = "none")+
  ggtitle("WT vs MUT: MetaCyc data") +
  theme(plot.title = element_text(hjust = 0.5))
baseline_pwy_mut

plot_grid(baseline_pwy_het, baseline_pwy_mut, labels = c("C", "D"))

