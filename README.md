# BasisGraphicalLasso

This package contains code to perform the basis graphical lasso
analysis of:

>Krock, M., Kleiber, W., and Becker, S. (2021), “Nonstationary modeling with sparsity for spatial data via the basis graphical lasso,” J. Comput. Graph. Statist., 30, 375–389, ISSN 1061-8600, URL https://doi.org/10.1080/10618600.2020.1811103

Please install the package with

```r
devtools::install_github("mlkrock/BasisGraphicalLasso")
```

The basic function is BGL, and the analysis of the minimum temperature
dataset in the paper can be reproduced as

```r
precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
     distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
     MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)
```
