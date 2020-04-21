# scry 1.0

Initial release of the scry package for analysis of high-dimensional data consisting of small counts (such as single-cell RNA-seq). Features included in this release:

* Feature selection with deviance
* Dimension reduction with GLM-PCA
* Approximate GLM-PCA using Pearson and deviance residuals from a null model.
* Optional adjustment for categorical batch effects for all methods.
* Support for SingleCellExperiment and SummarizedExperiment objects
* Support for matrix and sparse Matrix objects.
