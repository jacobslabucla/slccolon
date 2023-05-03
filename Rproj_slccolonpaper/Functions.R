## Make a taxa summary plot --
generate_L6_taxa_plots <- function(path_to_RDS, titlestring,greppattern, fillvector){
  #L2_lum<-readRDS("Long_Term/taxa_barplots/LuminalColon_level-6.RDS")
  #taxa <- gsub(".*g__","",taxa)
  #cols<-assign_cols
  titlestring<-c(titlestring)
  L2_lum<-readRDS(path_to_RDS)
  L2_lum<- as.matrix(L2_lum)
  L2_lum<-funrar::make_relative(L2_lum)
  L2_lum<-as.data.frame(t(L2_lum))
  toptaxa<- rowMeans(L2_lum)
  L2_lum$averageRA <-toptaxa
  L2_lum <- L2_lum %>% dplyr::mutate(keeptaxa = ifelse(averageRA >0.01, row.names(L2_lum), "Other"))
  L2_lum <-select(L2_lum,-averageRA)
  
  taxa<-L2_lum$keeptaxa
  L2_lum <- select(L2_lum,-keeptaxa)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  #taxa<-gsub(greppattern,"",taxa )
  
  family<-gsub(".*f__","",taxa )
  genus <- gsub(".*g__","",taxa)
  
  L2_lum$Family<-gsub(".g__.*","",family)
  L2_lum$Genus <- genus
  L2_lum <- L2_lum %>% mutate(annotation = ifelse(L2_lum$Genus=="", paste0(L2_lum$Family,"..f."), L2_lum$Genus))
  
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Family, Genus, annotation), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  L2_lum$Site <- factor(L2_lum$Site, levels=c("WT", "HET","MUT"))
  cols <- fillvector
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=annotation)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme(legend.position = "none")+
    theme_cowplot(20) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

## handle the genera names for taxonomy boxplot --
wrangle_genera_names <- function(csv_dataframe, filepathstring, rds_string){
  input_data<-read.csv(csv_dataframe, row.names=1, header=TRUE)
  
  taxa<-colnames(input_data)
  colnames <- strsplit(taxa, ".o__")
  
  order=new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(colnames)) {
    order[i] <- colnames[[i]][2]
    i=i+1
  }
  
  order<-unlist(order)
  order <- strsplit(order, ".f__")
  
  family =new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(order)) {
    family[i] <- order[[i]][2]
    i=i+1
  }
  
  order<-as.list(order)
  
  i=1
  for (i in 1:length(family)) {
    if (isFALSE(family[[i]]==".g__")) {
      family[[i]] = family[[i]] 
    }
    else {
      family[[i]] <- paste0(order[[i]]," (o)")   
      family[[i]] <- family[[i]][1]
    }
    i=i+1
  }
  
  
  family<-unlist(family)
  family <- strsplit(family, ".g__")
  
  genus =new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(family)) {
    genus[i] <- family[[i]][2]
    i=i+1
  }
  
  family<-as.list(family)
  
  i=1
  for (i in 1:length(genus)) {
    if (isFALSE(genus[[i]]=="NA")) {
      genus[[i]] = genus[[i]] 
    }
    else {
      
      genus[[i]] <- paste0(family[[i]]," (f)")   
    }
    i=i+1
  }
  
  

  colnames(input_data) <-as.character(genus)
  
  saveRDS(input_data, paste0(filepathstring, rds_string))
  
}

## Make a taxonomy dotplot --

make_taxa_dotplot <- function(ASV_significant_results_dataset, taxonomy_tsv_filepath, 
                              Relative_Abundance_filepath_rds,titlestring, colorvariable){
data<-as.data.frame(ASV_significant_results_dataset)
taxonomy <- read.delim(taxonomy_tsv_filepath)
taxonomy$feature <- taxonomy$Feature.ID
data <- merge(data,taxonomy, by="feature")
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub(";.*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub(";.*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub(";.*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)
data$Genus <-  gsub(";.*","",data$Genus)
data$Species <- gsub(".*s__","",data$Taxon)
data$annotation <- paste0(data$Genus," ", data$Species)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

#append relative abundance data 
relA <- readRDS(Relative_Abundance_filepath_rds)
relA$feature <- row.names(relA)
relA$Relative_Abundance <- relA$V1
data<-merge(data,relA,by="feature")
print(summary(data$Relative_Abundance))
max(data$Relative_Abundance)

#make graph
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
baseline_DAT <- ggplot(data, aes(x = coef, y = annotation, color = {{colorvariable}})) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Relative Abundance",range = c(0.5,8),
                        limits=c(sqrt(0.000001),sqrt(0.3)),
                        breaks=c(sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),
                        labels=c("0.0001","0.001","0.01","0.1")) + 
  #scale_color_manual(name="Phylum", values = phyla_colors)+
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme_cowplot(16) +
  ggtitle(titlestring) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right") 
baseline_DAT 

}

calculate_rsjensen <- function(data){
  relativedata<-sweep(data,2,colSums(data),"/")
  
  dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    
    for(i in 1:matrixColSize) {
      for(j in 1:matrixColSize) { 
        resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                               as.vector(inMatrix[,j]))
      }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
  }
  ## From: http://enterotype.embl.de/enterotypes.html
  
  data.dist=dist.JSD(relativedata)
  
  # to export distance matrix
  x=data.dist
  y=as.matrix(x)
  
  return(y) 
}

generate_pcoA_plots <- function(distance_matrix, counts, metadata, title, colorvariable,colorvector, wa_scores_filepath){
  # calculate mds
  lumcol.dist <- distance_matrix
  lumcol_counts <- counts 
  lumcol_meta <- metadata
  
  mds <- cmdscale(lumcol.dist, eig = TRUE, x.ret = TRUE)
  
  mds_values <- mds$points
  wa_scores <- wascores(mds_values, t(lumcol_counts))
  wa_scores <- data.frame(sample = rownames(wa_scores),
                          x = wa_scores[,1],
                          y = wa_scores[,2])
  
  # isolate taxa with strongest contribution to principal coordinate axes
  n_taxa <- 10
  wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
  wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
  wa_scores_final <- rbind(wa_scores_1, wa_scores_2)
  write.csv(wa_scores_final, {{wa_scores_filepath}})
  
  # calculate percentage of variation that each mds axis accounts for
  mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
  
  # plot
  mds_data <- data.frame(sample = rownames(mds_values),
                         x = mds_values[,1],
                         y = mds_values[,2])
  
  #merge phenotypic data 
  lumcol_meta$sample <- lumcol_meta$SampleID
  mds_meta <- merge(mds_data, lumcol_meta, by = "sample")
  
  
  p<- ggplot(mds_meta, aes(x, y, colour={{colorvariable}})) + 
    geom_point(size=3) + 
    labs(x = paste("PC1(", mds_var_per[1], "%)",sep=""),
         y = paste("PC2(", mds_var_per[2], "%)",sep="")) +
    scale_colour_manual(name="",values={{colorvector}}) +
    theme_cowplot(16)+
    theme(legend.position="top",legend.justification = "center") +
    labs(title= paste0({{title}})) 
  return(p)
}


assign_letter <- function(x) {
  first_num <- substr(x, start = 4, stop = 4)
  letter <- switch(first_num,
                   "1" = "Oxidoreductases",
                   "2" = "Transferases",
                   "3" = "Hydrolases",
                   "4" = "Lyases",
                   "5"="Isomerases",
                   "6"="Ligases","unknown")
  return(letter)
}

## Make a taxonomy dotplot for genus level --

make_genus_level_taxa_dotplot <- function(ASV_significant_results_dataset,
                                          Relative_Abundance_filepath_rds,titlestring){
  data <- as.data.frame(ASV_significant_results_dataset)
  data <- data %>% filter(qval <0.25)
  data <- data %>% filter(metadata=="Genotype")
  data$Taxon <- data$feature
  data$Phylum <- gsub(".*p__","",data$Taxon)
  data$Phylum <- gsub("\\..*","",data$Phylum)
  data$Family<- gsub(".*f__","",data$Taxon)
  data$Family <-  gsub("\\..*","",data$Family)
  data$Order<- gsub(".*o__","",data$Taxon)
  data$Order <-  gsub("\\..*","",data$Order)
  data$Genus<- gsub(".*g__","",data$Taxon)
  
  data$annotation <- gsub("\\.E","E",data$Genus)
  data$annotation <- gsub("\\.","_",data$annotation)
  data$annotation <- gsub("__","_",data$annotation)
  #data$Genus <- gsub("\\..*","",data$Genus)
  data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
  data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))
  
  #append relative abundance data 
  #relA <- readRDS("Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS")
  relA <- readRDS(Relative_Abundance_filepath_rds)
  relA$feature <- row.names(relA)
  relA$feature <- gsub(";",".",relA$feature)
  relA$feature <- gsub(" ",".",relA$feature)
  relA$feature <- gsub("-",".",relA$feature)
  relA$feature <- gsub("\\[",".",relA$feature)
  relA$feature <- gsub("\\]",".",relA$feature)
  relA$Relative_Abundance <- relA$V1
  data<-merge(data,relA,by="feature")
  print(data$feature)
  print(summary(data$Relative_Abundance))
  max(data$Relative_Abundance)
  
  #make graph
  y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  data$annotation= factor(as.character(data$annotation), levels = names(y))
  baseline_DAT <- ggplot(data, aes(x = coef, y = annotation ,color=Phylum)) + 
    geom_point(aes(color=Phylum,size = sqrt(Relative_Abundance))) +
    geom_point(aes(size=sqrt(Relative_Abundance)),shape = 1,colour = "black")+
    scale_size_continuous(name="Relative Abundance",range = c(0.5,8),
                          limits=c(sqrt(0.000001),sqrt(0.6)),
                          breaks=c(sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),
                          labels=c("0.0001","0.001","0.01","0.1")) + 
    #scale_color_manual(name="Phylum", values = phyla_colors)+
    geom_vline(xintercept = 0) + 
    xlab(label="Log2 Fold Change")+
    ylab(label=NULL)+
    theme_cowplot(16) +
    ggtitle(titlestring) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "right") 
  baseline_DAT 
  
}
