# _Trios_
## Luminal Colon
### Females
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Litter     4    1.5553 0.38883  7.5608 0.45328 9.999e-05 ***
Site       2    0.1981 0.09903  1.9256 0.05772    0.0102 *  
Genotype   2    0.2894 0.14468  2.8134 0.08433 9.999e-05 ***
Residuals 27    1.3885 0.05143         0.40467              
Total     35    3.4312                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Males
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.5161 0.51607  9.4114 0.10052 9.999e-05 ***
Litter          4    1.6111 0.40277  7.3450 0.31381 9.999e-05 ***
Site            2    0.3096 0.15478  2.8227 0.06030 9.999e-05 ***
Genotype        2    0.3392 0.16961  3.0931 0.06608 9.999e-05 ***
Residuals      43    2.3579 0.05484         0.45929              
Total          52    5.1338                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Mucosal Colon
### Females
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Litter     4    1.4212 0.35529  3.0831 0.31352 9.999e-05 ***
Site       2    0.3401 0.17004  1.4755 0.07502   0.05249 .  
Genotype   2    0.2365 0.11823  1.0259 0.05216   0.40136    
Residuals 22    2.5353 0.11524         0.55930              
Total     30    4.5330                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Males
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    0.6900 0.69004  9.7199 0.11178 9.999e-05 ***
Litter          4    1.4664 0.36661  5.1640 0.23754 9.999e-05 ***
Site            2    0.7928 0.39640  5.5837 0.12842 9.999e-05 ***
Genotype        2    0.3133 0.15667  2.2069 0.05076 0.0009999 ***
Residuals      41    2.9107 0.07099         0.47150              
Total          50    6.1733                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

-----
# _Baseline

### Females
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2   Pr(>F)   
Genotype   2   0.24576 0.122880  1.7566 0.08672 0.009699 **
Residuals 37   2.58830 0.069954         0.91328            
Total     39   2.83406                  1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Males
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Genotype   2    0.2712 0.13562  1.5762 0.06686 0.0356 *
Residuals 44    3.7858 0.08604         0.93314         
Total     46    4.0570                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```