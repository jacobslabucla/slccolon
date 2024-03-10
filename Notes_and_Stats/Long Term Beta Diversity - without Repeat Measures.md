### With HETs
Luminal colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.5510 0.55095  6.2363 0.05955 9.999e-05 ***
Site       2    0.4391 0.21957  2.4853 0.04747     2e-04 ***
Genotype   2    0.8406 0.42032  4.7577 0.09086 9.999e-05 ***
Residuals 84    7.4210 0.08835         0.80212              
Total     89    9.2517                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R

```
Mucosal Colon
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.4018 0.40177  3.2823 0.03048    0.0016 ** 
Site       2    1.8038 0.90188  7.3681 0.13685 9.999e-05 ***
Genotype   2    0.8159 0.40797  3.3330 0.06190    0.0003 ***
Residuals 83   10.1595 0.12240         0.77077              
Total     88   13.1810                 1.00000              
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
Sex        1    0.5488 0.54882  6.6099 0.13247 9.999e-05 ***
Site       2    0.2028 0.10138  1.2210 0.04894    0.2143    
Genotype   1    0.3194 0.31940  3.8467 0.07709 9.999e-05 ***
Residuals 37    3.0721 0.08303         0.74150              
Total     41    4.1431                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1    0.6115 0.61149  5.6933 0.11023 9.999e-05 ***
Site       2    0.7609 0.38045  3.5421 0.13716 9.999e-05 ***
Genotype   1    0.3086 0.30862  2.8734 0.05563  0.006899 ** 
Residuals 36    3.8666 0.10741         0.69699              
Total     40    5.5476                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```