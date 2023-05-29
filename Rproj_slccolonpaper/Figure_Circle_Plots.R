library(circlize)

here::i_am("Rproj_slccolonpaper/Figure_Circle_Plots.R")

### Mucosal ---
## Trios
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") 
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
data <- data %>% mutate(coef_abs = abs(coef))
write.csv(data, "Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/MUT_vs_WT_circleplot.csv")
plot <- data %>% 
  select(c("feature", "X2","coef_abs"))
plot <- data %>% 
  select(c("description", "X2","coef_abs"))


mat <- plot 

circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.6,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 

# Mucosal Colon MUT depleted 
data<-read.table("Long_Term/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") 
annotation <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_pathway_abundance/annotated_pathway.tsv", row.names=1)
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
data <- merge(data,df_new, by="feature")
data <- data %>% mutate(coef_abs = abs(coef))
write.csv(data, "Long_Term/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/MUT_vs_WT_circleplot.csv")

plot <- data %>% 
  select(c("feature", "X2","coef_abs"))
plot <- data %>% 
 select(c("description", "X2","coef_abs"))

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

### Luminal ---
## Trios
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") 
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
data <- data %>% mutate(coef_abs = abs(coef))
write.csv(data, "Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/MUT_vs_WT_circleplot.csv")

plot <- data %>% 
  select(c("feature", "X2","coef_abs"))
plot <- data %>% 
  select(c("description", "X2","coef_abs"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.6,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term
data<-read.table("Long_Term/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") 
annotation <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_pathway_abundance/annotated_pathway.tsv", row.names=1)
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
data <- data %>% mutate(coef_abs = abs(coef))
write.csv(data, "Long_Term/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/MUT_vs_WT_circleplot.csv")

plot <- data %>% 
  select(c("feature", "X2","coef_abs"))
plot <- data %>% 
  select(c("description", "X2","coef_abs"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.6,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Baseline
data<-read.table("Baseline/differential_Pathway/Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") 
annotation <- read.delim("Baseline/starting_files/picrust2_output_Baseline_ASV_table_Silva_v138_1.qza/export_pathway_abundance/annotated_pwy.tsv", row.names=1)
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
data <- data %>% mutate(coef_abs = abs(coef))
write.csv(data, "Baseline/differential_Pathway/Baseline_JAX_PWY_Maaslin2_Sex_Genotype/MUT_vs_WT_circleplot.csv")

plot <- data %>% 
  select(c("feature", "X2","coef_abs"))
plot <- data %>% 
  select(c("description", "X2","coef_abs"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.6,
              cex=1)
}, bg.border = NA) 
circos.clear()