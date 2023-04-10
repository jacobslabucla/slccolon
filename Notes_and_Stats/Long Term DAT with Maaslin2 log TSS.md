##### Model 8: ASV ~ Site + Sex + Genotype + (1/MouseID) + (1/Litter)
```R
Warning messages:
1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  unable to evaluate scaled gradient
2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
3: Model failed to converge with 1 negative eigenvalue: -8.4e-01 
4: Model failed to converge with 1 negative eigenvalue: -1.6e+00 
5: Model failed to converge with 1 negative eigenvalue: -8.2e-01
```

##### Model 1: ASV ~ Site + Sex + Genotype + (1/MouseID) + (1/Breeder)
```R
Model failed to converge with 1 negative eigenvalue: -2.7e+00 
```

##### Model 2: ASV ~ Site + Sex + Genotype + (1/MouseID) 
- all good, no probs