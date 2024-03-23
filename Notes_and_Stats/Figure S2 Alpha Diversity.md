## Trios 
### Luminal Colon
```R
Fixed effects:  shannon_entropy ~ Sequencing_Run + Sex + Site + Genotype 
                                Value  Std.Error DF   t-value p-value
(Intercept)                  5.452738 0.18719671 57 29.128388  0.0000
Sequencing_RunSLCMicrobiome  0.013938 0.23129595 25  0.060262  0.9524
SexMale                      0.130669 0.18898322 25  0.691433  0.4957
SiteDistal_Colon            -0.790264 0.09613095 57 -8.220702  0.0000
SiteProximal_Colon          -0.265951 0.09613095 57 -2.766549  0.0076
GenotypeHET                 -0.108608 0.20710671 25 -0.524407  0.6046
GenotypeMUT                 -0.112915 0.20659038 25 -0.546565  0.5895
```
```R
Fixed effects:  observed_features ~ Sequencing_Run + Sex + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                 230.14651 10.431537 57 22.062570  0.0000
Sequencing_RunSLCMicrobiome  91.98094 12.521239 25  7.345993  0.0000
SexMale                      -4.89760 10.234952 25 -0.478517  0.6364
SiteDistal_Colon            -34.85621  6.795298 57 -5.129460  0.0000
SiteProximal_Colon           -6.05621  6.795298 57 -0.891235  0.3765
GenotypeHET                 -20.41046 11.219317 25 -1.819225  0.0809
GenotypeMUT                  -9.03333 11.174311 25 -0.808402  0.4265
```

### Mucosal Colon
```R
Fixed effects:  shannon_entropy ~ Sequencing_Run + Sex + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                  4.841713 0.2774928 50 17.448070  0.0000
Sequencing_RunSLCMicrobiome  0.703704 0.3125923 25  2.251189  0.0334
SexMale                      0.261955 0.2647018 25  0.989622  0.3318
SiteDistal_Colon            -0.681204 0.2456118 50 -2.773497  0.0078
SiteProximal_Colon           0.516674 0.2342621 50  2.205540  0.0320
GenotypeHET                 -0.139480 0.2824031 25 -0.493904  0.6257
GenotypeMUT                 -0.303453 0.2928924 25 -1.036057  0.3101
```

```R
Fixed effects:  observed_features ~ Sequencing_Run + Sex + Site + Genotype 
                                Value Std.Error DF   t-value p-value
(Intercept)                 240.08330  15.03084 50 15.972710  0.0000
Sequencing_RunSLCMicrobiome 125.49955  17.18467 25  7.302994  0.0000
SexMale                     -17.69961  14.50498 25 -1.220244  0.2338
SiteDistal_Colon              1.30777  12.66464 50  0.103262  0.9182
SiteProximal_Colon            5.13285  12.06432 50  0.425457  0.6723
GenotypeHET                 -10.67698  15.50947 25 -0.688417  0.4975
GenotypeMUT                   0.58666  16.03459 25  0.036587  0.9711
```