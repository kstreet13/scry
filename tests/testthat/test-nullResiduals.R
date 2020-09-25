context("Test nullResiduals")
set.seed(1234)

test_that("nullResiduals works", {
    ncells <- 100
    u <- matrix(rpois(2000, 5), ncol=ncells)
    v <- log2(u + 1)
    require(Matrix)
    m <- Matrix(u)
    require(DelayedArray)
    d <- DelayedArray(u)
    require(SingleCellExperiment)
    se <- SummarizedExperiment(assays=list(logcounts=v,counts=u))
    sce <- SingleCellExperiment(assays=list(logcounts=v,counts=u,sparse_counts=m,delayed_counts=d))
    
    #check no error with different object types
    outU <- nullResiduals(u)
    expect_is(outU,"matrix")
    outD <- nullResiduals(d)
    expect_is(outD, 'DelayedArray')
    expect_warning({outM <- nullResiduals(m)},
                   'substantially larger than the input')
    expect_is(outM,"Matrix")
    outSE <- nullResiduals(se,assay="counts")
    expect_true("binomial_deviance_residuals" %in% assayNames(outSE))
    
    #check all object inputs give same output
    expect_equivalent(as.matrix(outM), outU)
    expect_equivalent(as.matrix(outD), outU)
    expect_equivalent(assay(outSE,"binomial_deviance_residuals"),outU)
    
    # check outputs are the same for in-memory and out-of-memory methods
    # given all possible arguments
    outU <- nullResiduals(u, fam = 'binomial', type = 'pearson')
    outD <- nullResiduals(d, fam = 'binomial', type = 'pearson')
    expect_equivalent(as.matrix(outD), outU)
    outU <- nullResiduals(u, fam = 'poisson', type = 'deviance')
    outD <- nullResiduals(d, fam = 'poisson', type = 'deviance')
    expect_equivalent(as.matrix(outD), outU)
    outU <- nullResiduals(u, fam = 'poisson', type = 'pearson')
    outD <- nullResiduals(d, fam = 'poisson', type = 'pearson')
    expect_equivalent(as.matrix(outD), outU)
    
    
    #check different input parameters don't throw errors
    outSCE <- nullResiduals(sce, assay="sparse_counts", fam = "binomial", type = "pearson")
    expect_true("binomial_pearson_residuals" %in% assayNames(outSCE))
    outSCE <- nullResiduals(sce, assay="counts", fam = "poisson")
    expect_true("poisson_deviance_residuals" %in% assayNames(outSCE))
    outSCE <- nullResiduals(sce, assay= 3, fam = "poisson", type = "pearson")
    expect_true("poisson_pearson_residuals" %in% assayNames(outSCE))
    sizeFactors(sce) <- compute_size_factors(assay(sce, 'counts'))
    outSCE <- nullResiduals(sce, assay=2, fam = "binomial")
    expect_true("binomial_deviance_residuals" %in% assayNames(outSCE))
    #verify it's a dense matrix
    expect_is(assay(outSCE,"binomial_deviance_residuals"),"matrix")
    
    #check handling of batch argument
    expect_error(nullResiduals(sce, batch = rep(1, ncol(sce))),
                 ') is not TRUE')
    outSCE <- nullResiduals(sce, batch = factor(rep(1:2, length.out = ncol(sce))))
    expect_true("binomial_deviance_residuals" %in% assayNames(outSCE))
    
    # check bad input types
    expect_error(nullResiduals(list(a = 1:5, b = letters[1:10])),
                 'currently not supported')
    expect_error(nullResiduals(1:20),
                 'currently not supported')
    # we may want to support data.frames later
    # expect_error(nullResiduals(data.frame(a=1:5,b=2:6)),
    #              'currently not supported')
    
    # hidden functions with redundancies
    expect_error(scry:::.binomial_deviance_residuals(u, rep(.5,13), colSums(u)),
                 ') is not TRUE')
    expect_error(scry:::.binomial_deviance_residuals(u, rep(.5,nrow(u)), 
                                                     rep(10,13)),
                 ') is not TRUE')
})
