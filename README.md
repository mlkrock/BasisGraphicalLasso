# BasisGraphicalLasso

This package contains code to perform the basis graphical lasso
analysis of Krock, Kleiber and Becker (2020+).

The basic function is BGL, and the analysis of the minimum temperature
dataset in the paper can be reproduced as

```r
precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
     distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
     MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)
```
