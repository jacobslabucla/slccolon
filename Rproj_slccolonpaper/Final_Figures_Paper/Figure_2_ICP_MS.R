library(here)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggpubr)
library(dplyr)
library(tidyr)
#install.packages("gridGraphics")

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_2_ICP_MS.R")

### Data Preprocessing ---
df<- as.data.frame(readr::read_csv(here("ICPMS/Analysis_ICP_MS.csv")))
row.names(df) <- df$SampleID
df <- df %>% select(-c("SampleID"))

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-"0"
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Batch <- paste0(df$Genotype, "_",df$Batch)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

icpms_summary <- df %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

icpms_summary %>% kable
# Generate a version of a df without outlier samples 
# create detect outlier function

detect_outlier <- function(x) {
  
  # calculate first quantile
  Quantile1 <- quantile(x, probs=.25,na.rm = TRUE)
  # calculate third quantile
  Quantile3 <- quantile(x, probs=.75,na.rm=TRUE)
  # calculate inter quartile range
  IQR = Quantile3-Quantile1
  # return true or false
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}

# create remove outlier function
remove_outlier <- function(dataframe,
                           columns=names(dataframe)) {
  
  # for loop to traverse in columns vector
  for (col in columns) {
    
    # remove observation if it satisfies outlier function
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  
  # return dataframe
  print("Remove outliers")
  print(dataframe)
  
}


# Subset by SampleType - without outliers
df_fp_col <- df %>% filter(SampleType=="FP-COL") 
df_fp_si <- df %>% filter(SampleType=="FP-SI") 
df_muc_col <- df %>% filter(SampleType=="MUC-COL") 
df_muc_si <- df %>% filter(SampleType=="MUC-SI") 
df_ts_col <- df %>% filter(SampleType=="TS-COL") 
df_ts_si <- df %>% filter(SampleType=="TS-SI") 

### Figures ---
generate_violin_plots <- function (input_data, X, titlestring) {
  # read in file
  data<-as.data.frame(input_data)
  
  #Ensure correct ordering of levels 
  data$Genotype <- factor(data$Genotype, levels = c("WT","MUT"))
  data$Genotype_Sex <- plyr::revalue(data$Genotype_Sex,c("WT_Male"="WT_M" ,
                                                         "WT_Female"="WT_F",
                                                         "MUT_Female"="MUT_F",
                                                         "MUT_Male"="MUT_M"))
  data$Genotype_Sex <- factor(data$Genotype_Sex, levels=c("WT_M","MUT_M", "WT_F", "MUT_F"))
  
  ggplot(data=data,aes(x={{X}},y=Concentration, fill=Genotype)) + 
    geom_boxplot(alpha=0.25,position=position_dodge(width=.75),size=1,color="black")+
    #scale_shape_manual(values=c(16,10))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25, height=0),alpha=0.8)+
    theme_cowplot(16) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(titlestring)+
    ylab("ug/g") +
    xlab("")
  #ylim(min,max) +
  
}


fp_col_plots <- list()
fp_si_plots <- list()
muc_col_plots <- list()
muc_si_plots <- list()
ts_col_plots <- list()
ts_si_plots <- list()


elementvector <- c("Iron", "Cobalt", "Copper", "Zinc", "Cadmium", "Manganese", "Selenium")

for (int in 1:7){
  print(int)
  element <- elementvector[int]
  fp_col <- df_fp_col %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element) 
  fp_si <- df_fp_si %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
  muc_col <- df_muc_col %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
  muc_si <- df_muc_si %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
  ts_col <- df_ts_col  %>%
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
  ts_si <- df_ts_si %>%
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
  
  fp_col_plots[[int]] <- fp_col
  fp_si_plots[[int]] <- fp_si
  muc_col_plots[[int]] <- muc_col
  muc_si_plots[[int]] <- muc_si
  ts_col_plots[[int]] <- ts_col
  ts_si_plots[[int]] <- ts_si
  
}

## Assemble Multi-Panel Figure
top_half <- plot_grid(fp_col_plots[[1]],fp_col_plots[[2]],
          fp_col_plots[[3]],fp_col_plots[[4]],
          fp_col_plots[[5]],fp_col_plots[[6]],
          nrow = 1, ncol=6,
          labels=c("A", "","","","",""
                   ),label_size = 20) 

cadmium <- muc_col_plots[[5]] + ylim(0,0.05)
middle <- plot_grid(muc_col_plots[[1]],muc_col_plots[[2]],
                      muc_col_plots[[3]],muc_col_plots[[4]],
                      cadmium,muc_col_plots[[6]],
                      nrow = 1, ncol=6,
                      labels=c(
                               "B","","","","",""),label_size = 20) 

#MUC COL
bottom_half <- plot_grid(ts_col_plots[[1]],ts_col_plots[[2]],
          ts_col_plots[[3]],ts_col_plots[[4]],
          ts_col_plots[[5]],ts_col_plots[[6]],
          nrow = 1, ncol=6, label_size = 20,
          labels=c("C", "","","","",""))



#MUC COL
middle_half <- plot_grid(muc_col_plots[[1]],muc_col_plots[[2]],
                    muc_col_plots[[3]],muc_col_plots[[4]],
                    cadmium,muc_col_plots[[6]],
                    nrow = 1, ncol=6,
                    labels=c(
                      "B","","","","",""),label_size = 20) 


#Plot arranged figures
plot_grid(plotlist = list(top_half, middle_half, bottom_half), nrow = 3)

### Add breakdown by genotype_sex variable ---

fp_col_plots <- list()
fp_si_plots <- list()
muc_col_plots <- list()
muc_si_plots <- list()
ts_col_plots <- list()
ts_si_plots <- list()


elementvector <- c("Iron", "Cobalt", "Copper", "Zinc", "Cadmium", "Manganese", "Selenium")

for (int in 1:7){
  print(int)
  element <- elementvector[int]
  fp_col <- df_fp_col %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  fp_si <- df_fp_si %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  muc_col <- df_muc_col %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  muc_si <- df_muc_si %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  ts_col <- df_ts_col  %>%
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  ts_si <- df_ts_si %>%
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  
  fp_col_plots[[int]] <- fp_col
  fp_si_plots[[int]] <- fp_si
  muc_col_plots[[int]] <- muc_col
  muc_si_plots[[int]] <- muc_si
  ts_col_plots[[int]] <- ts_col
  ts_si_plots[[int]] <- ts_si
  
}

## Assemble Multi-Panel Figure
top_half <- plot_grid(fp_col_plots[[1]],fp_col_plots[[2]],
                      fp_col_plots[[3]],fp_col_plots[[4]],
                      fp_col_plots[[5]],fp_col_plots[[6]],
                      nrow = 1, ncol=6,
                      labels=c("A", "","","","",""
                      ),label_size = 20) 
dev.new(width=10,height=10)
top_half

cadmium <- muc_col_plots[[5]] + ylim(0,0.05)
middle <- plot_grid(muc_col_plots[[1]],muc_col_plots[[2]],
                    muc_col_plots[[3]],muc_col_plots[[4]],
                    cadmium,muc_col_plots[[6]],
                    nrow = 1, ncol=6,
                    labels=c(
                      "B","","","","",""),label_size = 20) 
dev.new(width=10,height=10)
middle

bottom_half <- plot_grid(ts_col_plots[[1]],ts_col_plots[[2]],
                         ts_col_plots[[3]],ts_col_plots[[4]],
                         ts_col_plots[[5]],ts_col_plots[[6]],
                         nrow = 1, ncol=6, label_size = 20,
                         labels=c("C", "","","","",""))

dev.new(width=10,height=10)
bottom_half
### Statistics On Full Dataset--- 
df<- as.data.frame(readr::read_csv(here("ICPMS/Analysis_ICP_MS.csv")))
row.names(df) <- df$SampleID
df <- df %>% select(-c("SampleID"))

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-"0"
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Batch <- paste0(df$Genotype, "_",df$Batch)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

#Females 
df_fp_col <- df %>% filter(SampleType=="FP-COL" & Sex=="Female")
df_fp_si <- df %>% filter(SampleType=="FP-SI" & Sex=="Female")
df_muc_col <- df %>% filter(SampleType=="MUC-COL" & Sex=="Female")
df_muc_si <- df %>% filter(SampleType=="MUC-SI"& Sex=="Female")
df_ts_col <- df %>% filter(SampleType=="TS-COL"& Sex=="Female")
df_ts_si <- df %>% filter(SampleType=="TS-SI"& Sex=="Female")


element_stats_para <- list()
element_stats_nonpara <-list()


for (int in 1:7){
  print(int)
  ts_si_para <- t.test(df_ts_si[,int]~Genotype,df_ts_si)
  ts_si_nonpara <- wilcox.test(df_ts_si[,int]~Genotype,df_ts_si)
  ts_col_para <- t.test(df_ts_col[,int]~Genotype,df_ts_col)
  muc_si_para <- t.test(df_muc_si[,int]~Genotype,df_muc_si)
  muc_si_nonpara <- wilcox.test(df_muc_si[,int]~Genotype,df_muc_si)
  ts_col_nonpara <- wilcox.test(df_ts_col[,int]~Genotype,df_ts_col)
  muc_col_para <- t.test(df_muc_col[,int]~Genotype,df_muc_col)
  muc_col_nonpara <- wilcox.test(df_muc_col[,int]~Genotype,df_muc_col)
  fp_si_para <- t.test(df_fp_si[,int]~Genotype,df_fp_si)
  fp_si_nonpara <- wilcox.test(df_fp_si[,int]~Genotype,df_fp_si)
  fp_col_para <- t.test(df_fp_col[,int]~Genotype,df_fp_col)
  fp_col_nonpara <- wilcox.test(df_fp_col[,int]~Genotype,df_fp_col)
  
  element_stats_para[[int]] <-list(print(fp_col_para), print(fp_si_para),print(muc_col_para),print(muc_si_para), print(ts_col_para),print(ts_si_para))
  element_stats_nonpara[[int]] <-list(print(fp_col_nonpara), print(fp_si_nonpara),print(muc_col_nonpara),print(muc_si_nonpara), print(ts_col_nonpara),print(ts_si_nonpara))
  
}

# Iron 
element_stats_para[[1]]
element_stats_nonpara[[1]]

# Cobalt 
element_stats_para[[2]]
element_stats_nonpara[[2]]

# Copper
element_stats_para[[3]]
element_stats_nonpara[[3]]

# Zinc
element_stats_para[[4]]
element_stats_nonpara[[4]]

# Cadmium
element_stats_para[[5]]
element_stats_nonpara[[5]]

# Manganese
element_stats_para[[6]]
element_stats_nonpara[[6]]

# Selenium
element_stats_para[[7]]
element_stats_nonpara[[7]]

# Males
df_fp_col <- df %>% filter(SampleType=="FP-COL" & Sex=="Male")
df_fp_si <- df %>% filter(SampleType=="FP-SI" & Sex=="Male")
df_muc_col <- df %>% filter(SampleType=="MUC-COL" & Sex=="Male")
df_muc_si <- df %>% filter(SampleType=="MUC-SI"& Sex=="Male")
df_ts_col <- df %>% filter(SampleType=="TS-COL"& Sex=="Male")
df_ts_si <- df %>% filter(SampleType=="TS-SI"& Sex=="Male")


element_stats_para <- list()
element_stats_nonpara <-list()


for (int in 1:7){
  print(int)
  ts_si_para <- t.test(df_ts_si[,int]~Genotype,df_ts_si)
  ts_si_nonpara <- wilcox.test(df_ts_si[,int]~Genotype,df_ts_si)
  ts_col_para <- t.test(df_ts_col[,int]~Genotype,df_ts_col)
  muc_si_para <- t.test(df_muc_si[,int]~Genotype,df_muc_si)
  muc_si_nonpara <- wilcox.test(df_muc_si[,int]~Genotype,df_muc_si)
  ts_col_nonpara <- wilcox.test(df_ts_col[,int]~Genotype,df_ts_col)
  muc_col_para <- t.test(df_muc_col[,int]~Genotype,df_muc_col)
  muc_col_nonpara <- wilcox.test(df_muc_col[,int]~Genotype,df_muc_col)
  fp_si_para <- t.test(df_fp_si[,int]~Genotype,df_fp_si)
  fp_si_nonpara <- wilcox.test(df_fp_si[,int]~Genotype,df_fp_si)
  fp_col_para <- t.test(df_fp_col[,int]~Genotype,df_fp_col)
  fp_col_nonpara <- wilcox.test(df_fp_col[,int]~Genotype,df_fp_col)
  
  element_stats_para[[int]] <-list(print(fp_col_para), print(fp_si_para),print(muc_col_para),print(muc_si_para), print(ts_col_para),print(ts_si_para))
  element_stats_nonpara[[int]] <-list(print(fp_col_nonpara), print(fp_si_nonpara),print(muc_col_nonpara),print(muc_si_nonpara), print(ts_col_nonpara),print(ts_si_nonpara))
  
}

# Iron 
element_stats_para[[1]]
element_stats_nonpara[[1]]

# Cobalt 
element_stats_para[[2]]
element_stats_nonpara[[2]]

# Copper
element_stats_para[[3]]
element_stats_nonpara[[3]]

# Zinc
element_stats_para[[4]]
element_stats_nonpara[[4]]

# Cadmium
element_stats_para[[5]]
element_stats_nonpara[[5]]

# Manganese
element_stats_para[[6]]
element_stats_nonpara[[6]]

# Selenium
element_stats_para[[7]]
element_stats_nonpara[[7]]


