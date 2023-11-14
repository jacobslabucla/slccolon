library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(circlize)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_PWY_Maaslin2.R")

metadata <- read.table("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)
pathway <- read.delim("Trios/differential_Pathway/feature-table.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))
pathway <- pathway %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Grab samples in the pathway table that are present in counts --
names(pathway)
pathway <- pathway %>% select(c(names(counts)))

### Luminal SI ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Subset =="Luminal_SI", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum","Jejunum", "Ileum"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Trios/differential_Pathway/PICRUST2_PWY_Luminal_SI_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID"), 
                    fixed_effects = c("Sequencing_Run","Site",  "Sex", "Genotype"), normalization = "TSS", 
                    random_effects = c("MouseID"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

### Mucosal SI ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Subset =="Mucosal_SI", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum","Jejunum", "Ileum"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_SI_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID"), 
                    fixed_effects = c("Sequencing_Run","Site",  "Sex", "Genotype"), normalization = "TSS", 
                    random_effects = c("MouseID"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

### Visualization of the Results ---

create_bar_plot <- function(res_plot, title) {
  y <- tapply(res_plot$coef, res_plot$description, function(y) mean(y))
  y <- sort(y, FALSE)
  res_plot$description <- factor(as.character(res_plot$description), levels = names(y))
  cols <- c(WT = "black", HET = "blue", MUT = "firebrick")
  res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.25, abs(coef) > 0) %>%
    ggplot(aes(coef, description, fill = site)) +
    geom_bar(stat = "identity") +
    theme_cowplot(16) +
    theme(axis.text.y = element_text(face = "bold")) +
    scale_fill_manual(values = cols) +
    labs(x = "Effect size", y = "", fill = "") +
    theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

## Model 2: PWY ~ Sequencing_Run + Site + Sex + Genotype + (1|mouseID) --
# Mucosal SI - no significant differences by Genotype with q <0.25
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_SI_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")

res_plot <- data %>%
  filter(value == "MUT") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "MUT"))
pwy_mut <- create_bar_plot(res_plot, "Mucosal SI: WT vs MUT q<0.25")
res_plot <- data %>%
  filter(value == "HET") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Mucosal SI: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

# Luminal sI - no significant differences by Genotype with q <0.25 between WT and MUT but only between WT and HET 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_SI_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>%
  filter(value == "MUT") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "MUT"))
pwy_mut <- create_bar_plot(res_plot, "Luminal Colon: WT vs MUT q<0.25")
res_plot <- data %>%
  filter(value == "HET") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Luminal Colon: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

### Visualizing results as circle plots ---

# Mucosal SI MUT enriched 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
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
data <- data %>% 
  select(c("feature", "X1","coef"))

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

# Mucosal Colon MUT depleted 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
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
data <- data %>% 
  select(c("feature", "X1","coef"))
data <- unique(data)

mat <- data 
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