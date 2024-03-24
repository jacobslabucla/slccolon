# slccolon
This repository contains the original count data and analysis files used to produce the figures in the paper.

## Steps to Reproduce the Figures 
1. First, clone the repo. 
```bash
git clone https://github.com/jacobslabucla/slccolon
```
2. Then, open the `Rproj_slccolonpaper` folder. Go to `Final_Figures_Paper`.

3. Run everything in `Functions.R` first. 

4. Now, by running each Figure script (matches the Figures in the paper) you should be able to reproduce the figures in the paper.

## ICP MS Data 
- ICP-MS / Colon samples from 5-month old mice (WT, MUT) - n = 25, 23

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |MUT      |      12|
|Female |WT       |      13|
|Male   |MUT      |      11|
|Male   |WT       |      12|

- ICP-MS / Blood - n = 20,20 

|Sex    |Genotype | SampleID|
|:------|:--------|--------:|
|Female |MUT      |       10|
|Female |WT       |       10|
|Male   |MUT      |       10|
|Male   |WT       |       10|

## Microbiome Data
- SLC Long Term / 12-month old mice (WT, HET, MUT) - n= 7, 16, 7

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |HET      |      10|
|Female |MUT      |       6|
|Female |WT       |       4|
|Male   |HET      |       6|
|Male   |MUT      |       1|
|Male   |WT       |       3|

- SLC Trios / 2-month old mice (WT, HET, MUT) - n = 10, 10, 10
  
|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |HET      |       4|
|Female |MUT      |       4|
|Female |WT       |       4|
|Male   |HET      |       6|
|Male   |MUT      |       6|
|Male   |WT       |       6|

- Baseline / 3-4 month old mice (WT, HET, MUT) - n = 26, 30, 31

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |HET      |      10|
|Female |MUT      |      14|
|Female |WT       |      16|
|Male   |HET      |      20|
|Male   |MUT      |      17|
|Male   |WT       |      11|
  
## Histology Data 

- 5 months 

|Sex    |Genotype | Score|
|:------|:--------|-----:|
|Female |MUT      |     4|
|Female |WT       |     7|
|Male   |MUT      |     9|
|Male   |WT       |     9|


- 10 months 

|Sex    |Genotype | Score|
|:------|:--------|-----:|
|Female |HET      |     7|
|Female |MUT      |     5|
|Female |WT       |     4|
|Male   |HET      |     6|
|Male   |MUT      |    10|
|Male   |WT       |     4|

- Mucin 

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |HET      |       5|
|Female |MUT      |       7|
|Female |WT       |       6|
|Male   |HET      |       6|
|Male   |MUT      |       4|
|Male   |WT       |       6|

## Disease model Data 

- DSS 

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |WT       |       4|
|Female |HET      |       4|
|Female |MUT      |       3|
|Male   |WT       |       9|
|Male   |HET      |       9|
|Male   |MUT      |      12|

- TNBS

|Sex    |Genotype | MouseID|
|:------|:--------|-------:|
|Female |WT       |       8|
|Female |HET      |       8|
|Female |MUT      |       5|
|Male   |WT       |      12|
|Male   |HET      |      14|
|Male   |MUT      |      12|
