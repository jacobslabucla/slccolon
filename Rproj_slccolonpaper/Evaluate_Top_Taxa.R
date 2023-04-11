
## Trios Top Taxa --
data<-read.csv("Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv", header=TRUE,row.names=1)
data$feature <- data$sample
annotation <- read.delim("Trios/starting_files/final_taxonomy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","Taxon"))

data <- merge(data,annotation, by="feature")

data$Family <- gsub(".*; f__","",data$Taxon)
