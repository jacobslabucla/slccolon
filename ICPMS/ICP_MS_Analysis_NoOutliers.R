###### SLC ICP-MS Data Aggregated ---
### Date: 10.26.2022
### Contents: Includes data from 10/5/21 and 10/25/22 Core submissions
###### End ---

library(here)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/SLC Spontaneous/ICPMS/")
here::i_am("ICP_MS_Analysis.R")

### Data Preprocessing ---
df<- read.csv("Analysis_ICP_MS.csv", header=TRUE, row.names=1)

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-0
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Batch <- paste0(df$Genotype, "_",df$Batch)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

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
  data$Genotype_Batch <- factor(data$Genotype_Batch, levels=c("WT_One","WT_Two", "WT_Three", "MUT_One", "MUT_Two", "MUT_Three"))
  data$Genotype_Sex <- factor(data$Genotype_Sex, levels=c("WT_Male","WT_Female", "MUT_Male", "MUT_Female"))
  
  ggplot(data=data,aes(x={{X}},y=Concentration, fill=Genotype)) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    #scale_shape_manual(values=c(16,10))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=0.8)+
    theme_cowplot(16) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(titlestring)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
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

#FP COL
dev.new(width=15, height=10)
plot_grid(fp_col_plots[[1]],fp_col_plots[[2]],
          fp_col_plots[[3]],fp_col_plots[[4]],
          fp_col_plots[[5]],fp_col_plots[[6]],
          fp_col_plots[[7]],nrow = 2, ncol=4)

#FP SI
dev.new(width=15, height=10)
plot_grid(fp_si_plots[[1]],fp_si_plots[[2]],
          fp_si_plots[[3]],fp_si_plots[[4]],
          fp_si_plots[[5]],fp_si_plots[[6]],
          fp_si_plots[[7]],nrow = 2, ncol=4)

#MUC COL
dev.new(width=15, height=10)
plot_grid(muc_col_plots[[1]],muc_col_plots[[2]],
          muc_col_plots[[3]],muc_col_plots[[4]],
          muc_col_plots[[5]],muc_col_plots[[6]],
          muc_col_plots[[7]],nrow = 2, ncol=4)

#MUC SI
dev.new(width=15, height=10)
plot_grid(muc_si_plots[[1]],muc_si_plots[[2]],
          muc_si_plots[[3]],muc_si_plots[[4]],
          muc_si_plots[[5]],muc_si_plots[[6]],
          muc_si_plots[[7]],nrow = 2, ncol=4)

#TS COL
dev.new(width=15, height=10)
plot_grid(ts_col_plots[[1]],ts_col_plots[[2]],
          ts_col_plots[[3]],ts_col_plots[[4]],
          ts_col_plots[[5]],ts_col_plots[[6]],
          ts_col_plots[[7]],nrow = 2, ncol=4)

#TS SI
dev.new(width=15, height=10)
plot_grid(ts_si_plots[[1]],ts_si_plots[[2]],
          ts_si_plots[[3]],ts_si_plots[[4]],
          ts_si_plots[[5]],ts_si_plots[[6]],
          ts_si_plots[[7]],nrow = 2, ncol=4)