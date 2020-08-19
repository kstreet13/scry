context("Test nullResiduals")
set.seed(1234)

test_that("nullResiduals works", {
    ncells <- 100
    u <- matrix(rpois(2000, 5), ncol=ncells)
    v <- log2(u + 1)
    require(Matrix)
    m <- Matrix(u)
    require(SingleCellExperiment)
    se <- SummarizedExperiment(assays=list(logcounts=v,counts=u))
    sce <- SingleCellExperiment(assays=list(logcounts=v,counts=u,sparse_counts=m))
    
    #check no error with different object types
    outU<-nullResiduals(u)
    expect_is(outU,"matrix")
    outM<-nullResiduals(m)
    expect_is(outM,"matrix")
    outSE<-nullResiduals(se,assay="counts")
    expect_true("binomial_deviance_residuals" %in% assayNames(outSE))
    
    #check all object inputs give same output
    expect_equivalent(outM,outU)
    expect_equivalent(assay(outSE,"binomial_deviance_residuals"),outU)
    
    #check different input parameters don't throw errors
    outSCE <- nullResiduals(sce, assay="sparse_counts", fam = "binomial", type = "pearson")
    expect_true("binomial_pearson_residuals" %in% assayNames(outSCE))
    outSCE <- nullResiduals(sce, assay="counts", fam = "poisson")
    expect_true("poisson_deviance_residuals" %in% assayNames(outSCE))
    outSCE <- nullResiduals(sce, assay= 3, fam = "poisson", type = "pearson")
    expect_true("poisson_pearson_residuals" %in% assayNames(outSCE))
    outSCE <- nullResiduals(sce, assay=2, fam = "binomial")
    expect_true("binomial_deviance_residuals" %in% assayNames(outSCE))
    #verify it's a dense matrix
    expect_is(assay(outSCE,"binomial_deviance_residuals"),"matrix")
    
    #check SCE output equivalent to matrix output
    expect_equivalent(assay(outSCE,"binomial_deviance_residuals"),outU)
    
    #check handling of batch argument
    expect_error(nullResiduals(sce, batch = rep(1, ncol(sce))),
                 ') is not TRUE')
    outSCE <- nullResiduals(sce, batch = factor(rep(1:2, length.out = ncol(sce))))
    expect_true("binomial_deviance_residuals" %in% assayNames(outSCE))
})
