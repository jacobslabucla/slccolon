## Sex + Genotype
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Line        2    2.4757 1.23785 13.3004 0.16682 9.999e-05 ***
Sex         1    0.4153 0.41529  4.4622 0.02798 9.999e-05 ***
Genotype    2    0.2230 0.11151  1.1981 0.01503    0.2263    
Residuals 126   11.7267 0.09307         0.79017              
Total     131   14.8407                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
## Sex- genotype interaction 
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

              Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Line           2    2.4757 1.23785 13.3549 0.16682 9.999e-05 ***
Sex            1    0.4153 0.41529  4.4805 0.02798 9.999e-05 ***
Genotype       2    0.2230 0.11151  1.2031 0.01503    0.2220    
Sex:Genotype   2    0.2333 0.11663  1.2583 0.01572    0.1791    
Residuals    124   11.4934 0.09269         0.77445              
Total        131   14.8407                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```