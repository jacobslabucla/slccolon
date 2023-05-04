## Trios 
### With HETs (no repeat measures) - seed(11)
Luminal Colon 
```R

> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.6142 0.61419 10.2455 0.06531 9.999e-05 ***
Litter          7    3.2284 0.46119  7.6934 0.34330 9.999e-05 ***
Sex             1    0.3191 0.31909  5.3229 0.03393 9.999e-05 ***
Site            2    0.4691 0.23455  3.9127 0.04989 9.999e-05 ***
Genotype        2    0.2770 0.13851  2.3105 0.02946     2e-04 ***
Residuals      75    4.4960 0.05995         0.47811              
Total          88    9.4038                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.6142 0.61419  6.9644 0.06531 9.999e-05 ***
Sex             1    0.8231 0.82314  9.3337 0.08753 9.999e-05 ***
Site            2    0.4712 0.23558  2.6713 0.05010 9.999e-05 ***
Genotype        2    0.2637 0.13187  1.4953 0.02805    0.0457 *  
Residuals      82    7.2315 0.08819         0.76900              
Total          88    9.4038                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Mucosal Colon
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.9852 0.98517 10.2084 0.08393 9.999e-05 ***
Litter          7    2.6993 0.38562  3.9958 0.22996 9.999e-05 ***
Sex             1    0.3180 0.31801  3.2952 0.02709   0.00050 ***
Site            2    0.8953 0.44767  4.6388 0.07627 9.999e-05 ***
Genotype        2    0.2782 0.13912  1.4416 0.02370   0.05139 .  
Residuals      68    6.5624 0.09651         0.55905              
Total          81   11.7384                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.9852 0.98517  8.2016 0.08393 9.999e-05 ***
Sex             1    0.5055 0.50546  4.2080 0.04306 9.999e-05 ***
Site            2    0.9673 0.48367  4.0266 0.08241 9.999e-05 ***
Genotype        2    0.2716 0.13578  1.1304 0.02313    0.2496    
Residuals      75    9.0089 0.12012         0.76747              
Total          81   11.7384                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Without HETs 
Luminal Colon
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.4619 0.46192  8.4684 0.07265 9.999e-05 ***
Litter          7    2.6388 0.37698  6.9110 0.41501 9.999e-05 ***
Sex             1    0.2333 0.23334  4.2778 0.03670 9.999e-05 ***
Site            2    0.2914 0.14569  2.6710 0.04583 0.0004000 ***
Genotype        1    0.1692 0.16924  3.1026 0.02662 0.0006999 ***
Residuals      47    2.5637 0.05455         0.40320              
Total          59    6.3584                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.4619 0.46192  5.1390 0.07265 9.999e-05 ***
Sex             1    0.6086 0.60859  6.7707 0.09571 9.999e-05 ***
Site            2    0.2914 0.14569  1.6209 0.04583   0.03100 *  
Genotype        1    0.1427 0.14267  1.5873 0.02244   0.08069 .  
Residuals      54    4.8538 0.08989         0.76337              
Total          59    6.3584                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 11
```
Mucosal
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.6496 0.64964  6.7843 0.08860 9.999e-05 ***
Litter          7    1.9835 0.28336  2.9592 0.27051 9.999e-05 ***
Sex             1    0.3097 0.30973  3.2346 0.04224    0.0002 ***
Site            2    0.4822 0.24110  2.5179 0.06576 9.999e-05 ***
Genotype        1    0.1730 0.17299  1.8066 0.02359    0.0288 *  
Residuals      39    3.7345 0.09576         0.50930              
Total          51    7.3326                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.6496 0.64964  5.3845 0.08860 9.999e-05 ***
Sex             1    0.4143 0.41430  3.4339 0.05650 9.999e-05 ***
Site            2    0.5706 0.28531  2.3648 0.07782    0.0003 ***
Genotype        1    0.1481 0.14810  1.2276 0.02020    0.2060    
Residuals      46    5.5499 0.12065         0.75688              
Total          51    7.3326                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```