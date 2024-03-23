# _Trios_
## Luminal Colon

### Females
```R
Fixed effects:  shannon_entropy ~ Litter + Site + Genotype 
                Value  Std.Error DF  t-value p-value
(Intercept)  5.016097 0.11665463 22 42.99955  0.0000
Litter1     -0.307506 0.20058226  5 -1.53307  0.1858
Litter2      0.620781 0.20058226  5  3.09489  0.0270
Litter3     -0.312039 0.20058226  5 -1.55567  0.1805
Litter4      0.126278 0.25278420  5  0.49955  0.6386
Site1        0.274781 0.09187013 22  2.99097  0.0067
Site2       -0.345877 0.09187013 22 -3.76485  0.0011
Genotype1    0.023712 0.15503678  5  0.15294  0.8844
Genotype2   -0.122626 0.15503678  5 -0.79095  0.4648
```

```R
Fixed effects:  shannon_entropy ~ Site + Genotype 
                       Value Std.Error DF   t-value p-value
(Intercept)         5.346468 0.2757518 22 19.388697  0.0000
SiteDistal_Colon   -0.620658 0.1591236 22 -3.900480  0.0008
SiteProximal_Colon -0.203685 0.1591236 22 -1.280043  0.2139
GenotypeHET        -0.146338 0.3676926  9 -0.397989  0.6999
GenotypeMUT         0.011755 0.3676926  9  0.031969  0.9752
```

```R
Fixed effects:  observed_features ~ Litter + Site + Genotype 
                Value Std.Error DF  t-value p-value
(Intercept) 204.50370  5.453553 22 37.49917  0.0000
Litter1     -16.72593  9.377133  5 -1.78369  0.1346
Litter2      43.05185  9.377133  5  4.59115  0.0059
Litter3      -5.39259  9.377133  5 -0.57508  0.5901
Litter4     -15.57778 11.817551  5 -1.31819  0.2446
Site1         6.72222  5.535000 22  1.21449  0.2374
Site2       -15.69444  5.535000 22 -2.83549  0.0096
Genotype1     1.65741  7.247902  5  0.22867  0.8282
Genotype2    -6.84259  7.247902  5 -0.94408  0.3885
```

```R
Fixed effects:  observed_features ~ Site + Genotype 
                       Value Std.Error DF   t-value p-value
(Intercept)        214.22222 16.230629 22 13.198640  0.0000
SiteDistal_Colon   -22.41667  9.586887 22 -2.338263  0.0289
SiteProximal_Colon   2.25000  9.586887 22  0.234696  0.8166
GenotypeHET         -8.50000 21.577635  9 -0.393926  0.7028
GenotypeMUT          6.08333 21.577635  9  0.281928  0.7844
```
### Males 
```R
Fixed effects:  shannon_entropy ~ Litter + Site + Genotype 
                Value  Std.Error DF  t-value p-value
(Intercept)  5.165252 0.10539682 33 49.00767  0.0000
Litter1      0.298871 0.23513552 10  1.27106  0.2325
Litter2      0.256743 0.23835085 10  1.07716  0.3067
Litter3     -0.286510 0.23513552 10 -1.21849  0.2510
Litter4     -0.154769 0.23513552 10 -0.65821  0.5253
Litter5     -0.113709 0.23513552 10 -0.48359  0.6391
Site1        0.409821 0.06949871 33  5.89682  0.0000
Site2       -0.502847 0.06811585 33 -7.38223  0.0000
Genotype1    0.090050 0.14884067 10  0.60501  0.5587
Genotype2    0.015929 0.14947858 10  0.10656  0.9172
```
```R
Fixed effects:  shannon_entropy ~ Sequencing_Run + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                  5.658412 0.2096941 33 26.984118  0.0000
Sequencing_RunSLCMicrobiome  0.011444 0.2242720 14  0.051028  0.9600
SiteDistal_Colon            -0.908322 0.1194772 33 -7.602477  0.0000
SiteProximal_Colon          -0.312450 0.1194772 33 -2.615143  0.0133
GenotypeHET                 -0.078467 0.2596824 14 -0.302165  0.7670
GenotypeMUT                 -0.196028 0.2586086 14 -0.758012  0.4610
```
```R
Fixed effects:  observed_features ~ Litter + Site + Genotype 
                Value Std.Error DF  t-value p-value
(Intercept) 232.50426  4.303281 33 54.02954  0.0000
Litter1      79.27351  9.552652 10  8.29859  0.0000
Litter2     -37.03424  9.963981 10 -3.71681  0.0040
Litter3      43.27351  9.552652 10  4.53000  0.0011
Litter4     -28.28204  9.552652 10 -2.96065  0.0143
Litter5     -41.28204  9.552652 10 -4.32153  0.0015
Site1        18.34186  5.404661 33  3.39371  0.0018
Site2       -24.94871  5.311016 33 -4.69754  0.0000
Genotype1    15.77351  6.058214 10  2.60366  0.0263
Genotype2   -12.43592  6.140475 10 -2.02524  0.0704
```
```R
Fixed effects:  observed_features ~ Sequencing_Run + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                 236.18551 10.870638 33 21.726923  0.0000
Sequencing_RunSLCMicrobiome  91.78869 10.684682 14  8.590681  0.0000
SiteDistal_Colon            -43.53373  9.308075 33 -4.676985  0.0000
SiteProximal_Colon          -11.97817  9.308075 33 -1.286858  0.2071
GenotypeHET                 -27.96627 12.418024 14 -2.252071  0.0409
GenotypeMUT                 -19.11111 12.297202 14 -1.554102  0.1425
```
## Mucosal Colon

### Females
```R
Fixed effects:  shannon_entropy ~ Litter + Site + Genotype 
                Value Std.Error DF   t-value p-value
(Intercept)  4.740356 0.2988451 17 15.862253  0.0000
Litter1      0.203050 0.5509100  5  0.368573  0.7275
Litter2     -0.185184 0.4837070  5 -0.382844  0.7176
Litter3     -0.702131 0.5039235  5 -1.393329  0.2223
Litter4      0.271531 0.6058741  5  0.448163  0.6728
Site1        0.273782 0.3038979 17  0.900900  0.3802
Site2       -0.756197 0.3297112 17 -2.293515  0.0348
Genotype1    0.089831 0.3861087  5  0.232658  0.8253
Genotype2   -0.051301 0.3806723  5 -0.134764  0.8981
```
```R
Fixed effects:  shannon_entropy ~ Site + Genotype 
                       Value Std.Error DF   t-value p-value
(Intercept)         4.991247 0.4903890 17 10.178139  0.0000
SiteDistal_Colon   -1.131238 0.5331015 17 -2.121993  0.0488
SiteProximal_Colon  0.244930 0.5174904 17  0.473304  0.6420
GenotypeHET        -0.109940 0.5561675  9 -0.197674  0.8477
GenotypeMUT        -0.095639 0.6130681  9 -0.156001  0.8795
```
```R
Fixed effects:  observed_features ~ Litter + Site + Genotype 
                Value Std.Error DF   t-value p-value
(Intercept) 240.69271  12.43531 17 19.355586  0.0000
Litter1      -4.88286  23.11767  5 -0.211218  0.8411
Litter2      19.64062  19.66400  5  0.998811  0.3637
Litter3     -21.05873  20.76544  5 -1.014124  0.3571
Litter4     -16.37567  24.54587  5 -0.667145  0.5342
Site1         8.55543  14.80718 17  0.577789  0.5710
Site2       -16.99314  16.06457 17 -1.057802  0.3049
Genotype1    -1.86053  15.78538  5 -0.117864  0.9108
Genotype2    -8.77355  15.51997  5 -0.565307  0.5963
```
```R
Fixed effects:  observed_features ~ Site + Genotype 
                       Value Std.Error DF   t-value p-value
(Intercept)        240.45329  21.95189 17 10.953650  0.0000
SiteDistal_Colon   -26.33407  25.39567 17 -1.036951  0.3143
SiteProximal_Colon   3.75400  24.69817 17  0.151995  0.8810
GenotypeHET         -6.67660  23.91961  9 -0.279127  0.7865
GenotypeMUT         25.81673  26.76335  9  0.964630  0.3599
```
### Males
```R
Fixed effects:  shannon_entropy ~ Litter + Site + Genotype 
                Value Std.Error DF  t-value p-value
(Intercept)  5.127128 0.1363382 31 37.60596  0.0000
Litter1      0.650983 0.3001484 10  2.16887  0.0553
Litter2     -0.383439 0.3001484 10 -1.27750  0.2303
Litter3      0.304530 0.3001484 10  1.01460  0.3342
Litter4      0.130439 0.3001484 10  0.43458  0.6731
Litter5     -0.124610 0.3076642 10 -0.40502  0.6940
Site1       -0.088972 0.1272459 31 -0.69921  0.4896
Site2       -0.530979 0.1353887 31 -3.92189  0.0005
Genotype1    0.196450 0.1938818 10  1.01324  0.3348
Genotype2    0.057797 0.1909538 10  0.30267  0.7683
```
```R
Fixed effects:  shannon_entropy ~ Sequencing_Run + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                  5.003386 0.2744844 31 18.228308  0.0000
Sequencing_RunSLCMicrobiome  0.703228 0.2737052 14  2.569291  0.0223
SiteDistal_Colon            -0.415194 0.2296688 31 -1.807796  0.0804
SiteProximal_Colon           0.708922 0.2163604 31  3.276581  0.0026
GenotypeHET                 -0.150780 0.3189067 14 -0.472803  0.6436
GenotypeMUT                 -0.448137 0.3209567 14 -1.396254  0.1844
```
```R
Fixed effects:  observed_features ~ Litter + Site + Genotype 
                Value Std.Error DF   t-value p-value
(Intercept) 264.39724  8.645490 31 30.582099  0.0000
Litter1      87.49165 19.118062 10  4.576387  0.0010
Litter2     -67.39724 19.118062 10 -3.525317  0.0055
Litter3      77.04721 19.118062 10  4.030074  0.0024
Litter4     -34.06391 19.118062 10 -1.781766  0.1051
Litter5     -52.43993 19.463410 10 -2.694283  0.0225
Site1       -10.28613  6.753620 31 -1.523054  0.1379
Site2        11.07226  7.195696 31  1.538733  0.1340
Genotype1     9.58453 12.277336 10  0.780669  0.4531
Genotype2    -5.17502 12.142218 10 -0.426200  0.6790
```
```R
Fixed effects:  observed_features ~ Sequencing_Run + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                 222.21292  17.56769 31 12.648958  0.0000
Sequencing_RunSLCMicrobiome 124.06254  18.31714 14  6.773029  0.0000
SiteDistal_Colon             20.04159  12.15005 31  1.649507  0.1091
SiteProximal_Colon            9.50000  11.40563 31  0.832922  0.4113
GenotypeHET                 -14.19207  21.27586 14 -0.667050  0.5156
GenotypeMUT                 -14.17589  21.36691 14 -0.663451  0.5178
```
# _Baseline_
### Females
```R
Call:
lm(formula = observed_otus ~ Genotype, data = f_data_meta)

Residuals:
    Min      1Q  Median      3Q     Max 
-35.625  -9.343  -2.086  10.513  47.375 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  134.599      3.141  42.847   <2e-16 ***
Genotype1     -5.974      4.216  -1.417    0.165    
Genotype2      7.501      4.746   1.581    0.122   
```

```R
Call:
lm(formula = shannon ~ Genotype, data = f_data_meta)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.96489 -0.23055 -0.01355  0.29069  0.94106 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 4.977439   0.071864  69.262   <2e-16 ***
Genotype1   0.002636   0.096457   0.027    0.978    
Genotype2   0.052472   0.108571   0.483    0.632  
```

### Males
```R
Call:
lm(formula = shannon ~ Genotype, data = m_data_meta)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.86484 -0.21279  0.04336  0.30117  0.99062 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  5.05462    0.08826  57.270   <2e-16 ***
Genotype1    0.07526    0.13503   0.557    0.580    
Genotype2   -0.02827    0.11634  -0.243    0.809
```

```R
Call:
lm(formula = observed_otus ~ Genotype, data = m_data_meta)

Residuals:
    Min      1Q  Median      3Q     Max 
-85.250  -7.875  -0.455  14.750  73.500 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 138.40152    4.37872  31.608   <2e-16 ***
Genotype1     0.05303    6.69928   0.008    0.994    
Genotype2    -2.90152    5.77167  -0.503    0.618      
```