```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Line        2    99.452  49.726  39.256 0.37334 9.999e-05 ***
Sex         1     5.772   5.772   4.557 0.02167    0.0111 *  
Genotype    2     1.555   0.778   0.614 0.00584    0.6557    
Residuals 126   159.604   1.267         0.59915              
Total     131   266.384                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

              Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Line           2    99.452  49.726  38.777 0.37334 9.999e-05 ***
Sex            1     5.772   5.772   4.501 0.02167    0.0119 *  
Genotype       2     1.555   0.778   0.606 0.00584    0.6640    
Sex:Genotype   2     0.593   0.297   0.231 0.00223    0.9240    
Residuals    124   159.011   1.282         0.59693              
Total        131   266.384                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```