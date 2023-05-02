

here::i_am("Rproj_slccolonpaper/Figure_Circle_Enzymes_plots.R")

### Luminal ---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Baseline
data<-read.table("Baseline/differential_EC/Baseline_JAX_EC_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Mucosal ---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Luminal Depleted---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Baseline
data<-read.table("Baseline/differential_EC/Baseline_JAX_EC_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Mucosal ---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Mucosal Depleted---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Luminal all---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT")
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") 
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Baseline
data<-read.table("Baseline/differential_EC/Baseline_JAX_EC_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT")
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

### Mucosal all---
## Trios 
data<-read.table("Trios/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT") 
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

## Long Term 
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)

data <- data %>% 
  filter(value=="MUT")
data$enzyme_class <- sapply(data$feature, assign_letter)

plot <- data %>% 
  select(c("feature", "enzyme_class","coef"))

mat <- plot
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()
