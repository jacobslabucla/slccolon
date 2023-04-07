library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/SLC_GitHub/slcproject/SLC_Microbiome_Baseline/")

### Note: First remove "#Constructed from biom file row"
### Run Maaslin2 and get table of relative abundances

run_Maaslin2 <- function(counts_filepath, metadata_filepath, subset_string) {
  #input_data <- read.delim("export_s3_min10000_PFF_Baseline_min10000_no_tax_PFF_ASV_table/feature-table.tsv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
  input_data <- read.delim(counts_filepath, header=TRUE, row.names=1) # choose filtered non rarefied csv file
  
  
  df_input_data <- as.data.frame(input_data)
  df_input_data <- select(df_input_data, -c("taxonomy"))
  
  transposed_input_data <- t(df_input_data)
  transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
  df_relative_ASV <- make_relative(transposed_input_data)
  df_relative_ASV <- as.data.frame(df_relative_ASV)
  Relative_Abundance <- summarize_all(df_relative_ASV, mean)
  Relative_Abundance <- as.data.frame(t(Relative_Abundance))
  
  readr::write_rds(Relative_Abundance,paste0("Relative_Abundance_",subset_string,"_ASV.RDS"))
  
  input_metadata <-read.delim(metadata_filepath,sep="\t",header=TRUE, row.names=1)
  #input_metadata <-read.delim("../starting_files/PFF_Mapping.tsv",sep="\t",header=TRUE, row.names=1)
  
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
                      output = paste0("ASV-level",subset_string,"_Maaslin2_Line_Sex_Genotype"), 
                      fixed_effects = c("Line", "Sex","Genotype"),normalization="TSS", 
                      min_prevalence = 0.14,
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}


## Sex and Genotype

# Baseline
run_Maaslin2("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_ec_metagenome/feature-table.tsv",
             "starting_files/Baseline_Metadata - Baseline_Metadata.tsv","SLC_Baseline_EC")
run_Maaslin2("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_ko_metagenome/feature-table.tsv",
             "starting_files/Baseline_Metadata - Baseline_Metadata.tsv","SLC_Baseline_KO")
run_Maaslin2("starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_pathway_abundance/feature-table.tsv",
             "starting_files/Baseline_Metadata - Baseline_Metadata.tsv","SLC_Baseline_PWY")


### Visualize EC results  ---

# Baseline WT vs HET
data<-read.table("ASV-level_SLC_Baseline_EC_Maaslin2_Line_Sex_Genotype/significant_results.tsv", header=TRUE)
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

