library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(circlize)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Baseline_JAX_PWY_Maaslin2.R")
setwd("..")
metadata <- read.table("Baseline/starting_files/Baseline_Metadata.tsv", header=TRUE)
counts <- read.table("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)
pathway <- read.delim("Baseline/differential_Pathway/feature-table.tsv", header = TRUE, row.names=1)
enzyme <- read.delim("Baseline/differential_EC/feature-table.tsv", header = TRUE, row.names=1)
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

### Visualize KO results  ---

# Baseline WT vs HET
data<-read.table("Baseline/Baseline_JAX_KO_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Baseline/differential_KO/annotated_KO.tsv", row.names=1)
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
  filter(qval < 0.25, abs(coef) > 0) %>%
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
data<-read.table("Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("differential_Pathway/annotated_pwy.tsv", row.names=1)
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
data<-read.table("Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("differential_Pathway/annotated_pwy.tsv", row.names=1)
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
  filter(qval < 0.25, abs(coef) > 0) %>%
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


### Visualizing results as circle plots ---

# Baseline MUT enriched 
data<-read.table("Baseline/Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
annotation <- read.delim("Baseline/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

#Classification strategy no 1 
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0) %>%
  select(c("feature", "classification","coef"))

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
plot <- data %>% 
  select(c("description", "X2","coef"))

mat <- plot 
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.5,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Baseline MUT depleted 
data<-read.table("Baseline/Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
annotation <- read.delim("Baseline/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# Classification strategy no 1
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  select(c("feature", "classification","coef"))
data<- remove_missing(data)

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
plot <- data %>% 
  select(c("feature", "X1","coef"))

mat <- plot
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Luminal Colon MUT enriched 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)

# Merge higher_classification column
data <- merge(data,df_new, by="feature")
data <- data %>% 
  select(c("description", "X1","coef"))

mat <- data 
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.25,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Luminal Colon MUT depleted 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# Classification strategy no 1
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  select(c("feature", "classification","coef"))
data<- remove_missing(data)

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
plot <- data %>% 
  select(c("feature", "X1","coef"))
plot <- data %>% 
  select(c("feature", "X2","coef"))

mat <- plot
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()
