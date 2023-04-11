## JAX mice only 
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.3766 0.37656  4.6215 0.05081 9.999e-05 ***
Genotype   2    0.2714 0.13572  1.6657 0.03663  0.008999 ** 
Residuals 83    6.7628 0.08148         0.91256              
Total     86    7.4108                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Full background
```R
> data.adonis$aov.tab 
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Background   1    2.0513 2.05133 24.2513 0.15307 9.999e-05 ***
Sex          1    0.3632 0.36317  4.2935 0.02710 9.999e-05 ***
Genotype     2    0.2443 0.12213  1.4439 0.01823   0.05339 .  
Residuals  127   10.7425 0.08459         0.80160              
Total      131   13.4012                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 
```
