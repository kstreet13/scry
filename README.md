# scry

<!-- badges: start -->
[![R-CMD-check](https://github.com/kstreet13/scry/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/kstreet13/scry/actions)
[![codecov](https://codecov.io/gh/kstreet13/scry/branch/master/graph/badge.svg?token=2QCzltvkbJ)](https://codecov.io/gh/kstreet13/scry)
<!-- badges: end -->

<img src=inst/scry_sticker.png height="200">

The released version of scry can be found on Bioconductor at https://bioconductor.org/packages/scry.

A collection of methods for the analysis of small-count data, such as single-cell RNA-seq.

Included methods:
 - Deviance for feature selection
 - GLM-PCA for dimension reduction
 - Null residuals, a transformation that when combined with PCA approximates GLM-PCA

Planned additions for upcoming release:
 - Quasi-UMI normalization for single-cell read counts
 - Cell type prediction for scRNA-seq

## Installation
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("scry")
```

The Bioconductor release requires the latest version of R. If you prefer to use the older
R 3.6, you can install from the unofficial r3 branch as follows:

```
remotes::install_github("kstreet13/scry@r3")
```

The r3 branch is not tested as regularly.


## Issues and bug reports
Please use https://github.com/kstreet13/scry/issues to submit issues, bug reports, and comments.
