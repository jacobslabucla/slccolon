### With HET
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.3709 0.37090  4.6895 0.05152 9.999e-05 ***
Genotype   2    0.2638 0.13188  1.6674 0.03664  0.009699 ** 
Residuals 83    6.5646 0.07909         0.91184              
Total     86    7.1993                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Without HET
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.3373 0.33726  4.6573 0.07666 9.999e-05 ***
Genotype   1    0.1515 0.15149  2.0920 0.03444  0.005899 ** 
Residuals 54    3.9104 0.07241         0.88890              
Total     56    4.3992                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```