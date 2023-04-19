
## Trios Top Taxa --
data<-read.csv("Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv", header=TRUE,row.names=1)
data$feature <- data$sample
annotation <- read.delim("Trios/starting_files/final_taxonomy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","Taxon"))

data <- merge(data,annotation, by="feature")
data$Phylum <- gsub(".*; p__","",data$Taxon)
data$Phylum <- gsub(";.*","",data$Phylum)
data$Order<- gsub(".*; o__","",data$Taxon)
data$Order <- gsub(";.*","",data$Order)
data$Family<- gsub(".*; f__","",data$Taxon)
data$Family <- gsub(";.*","",data$Family)
data$Genus<- gsub(".*; g__","",data$Taxon)
data$Genus <- gsub("; s__"," ",data$Genus)

#data <- data %>% mutate(Taxon2 = ifelse(data$Species=="", paste(data$Genus,"(g)"), data$Species))
data <- data %>% mutate(label = ifelse(data$Genus==" ", paste(data$Family,"(f)"), data$Genus))
data <- data %>% mutate(Taxon = ifelse(data$label==" (f)", paste(data$Order,"(o)"), data$label))

ggplot(data, aes(x = x, y = Taxon)) +        # Create barchart with ggplot2
  geom_bar(stat = "identity") +
  theme_cowplot(16)+
  ggtitle("Trios: Luminal Colon") +
  ylab("Top ASV contributions to PC1")


## Long Term Top Taxa --
data<-read.csv("Long_Term/beta_diversity/LumCol_Top_Taxa_PcoA.csv", header=TRUE,row.names=1)
data$feature <- data$sample
annotation <- read.delim("Long_Term/starting_files/final_taxonomy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","Taxon"))

data <- merge(data,annotation, by="feature")
data$Phylum <- gsub(".*; p__","",data$Taxon)
data$Phylum <- gsub(";.*","",data$Phylum)
data$Order<- gsub(".*; o__","",data$Taxon)
data$Order <- gsub(";.*","",data$Order)
data$Family<- gsub(".*; f__","",data$Taxon)
data$Family <- gsub(";.*","",data$Family)
data$Genus<- gsub(".*; g__","",data$Taxon)
data$Genus <- gsub("; s__"," ",data$Genus)

#data <- data %>% mutate(Taxon2 = ifelse(data$Species=="", paste(data$Genus,"(g)"), data$Species))
data <- data %>% mutate(label = ifelse(data$Genus==" ", paste(data$Family,"(f)"), data$Genus))
data <- data %>% mutate(Taxon = ifelse(data$label==" (f)", paste(data$Order,"(o)"), data$label))

ggplot(data, aes(x = x, y = Taxon)) +        # Create barchart with ggplot2
  geom_bar(stat = "identity") +
  theme_cowplot(16)+
  ggtitle("Long Term: Luminal Colon") +
  ylab("Top ASV contributions to PC1")

## Baseline Top Taxa --
# JAX only Dataset
data<-read.csv("Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv", header=TRUE,row.names=1)
data$feature <- data$sample
annotation <- read.delim("Baseline/starting_files/Taxonomy_Key.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","Taxon"))

data <- merge(data,annotation, by="feature")
data$Phylum <- gsub(".*; p__","",data$Taxon)
data$Phylum <- gsub(";.*","",data$Phylum)
data$Order<- gsub(".*; o__","",data$Taxon)
data$Order <- gsub(";.*","",data$Order)
data$Family<- gsub(".*; f__","",data$Taxon)
data$Family <- gsub(";.*","",data$Family)
data$Genus<- gsub(".*; g__","",data$Taxon)
data$Genus <- gsub("; s__"," ",data$Genus)

#data <- data %>% mutate(Taxon2 = ifelse(data$Species=="", paste(data$Genus,"(g)"), data$Species))
data <- data %>% mutate(label = ifelse(data$Genus==" ", paste(data$Family,"(f)"), data$Genus))
data <- data %>% mutate(Taxon = ifelse(data$label==" (f)", paste(data$Order,"(o)"), data$label))

ggplot(data, aes(x = x, y = Taxon)) +        # Create barchart with ggplot2
  geom_bar(stat = "identity") +
  theme_cowplot(16)+
  ggtitle("Baseline JAX: Luminal Colon") +
  ylab("Top ASV contributions to PC1")

# Full Dataset
data<-read.csv("Baseline/beta_diversity/Baseline_Top_Taxa_PcoA.csv", header=TRUE,row.names=1)
data$feature <- data$sample
annotation <- read.delim("Baseline/starting_files/Taxonomy_Key.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","Taxon"))

data <- merge(data,annotation, by="feature")

data$Family <- gsub(".*; f__","",data$Taxon)
