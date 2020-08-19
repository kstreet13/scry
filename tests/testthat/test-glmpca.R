context("Test GLMPCA")
set.seed(1234)

test_that("GLMPCA works", {
    ncells <- 100
    u <- matrix(rpois(2000, 5), ncol=ncells)
    v <- log2(u + 1)
    
    require(SingleCellExperiment)
    se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
    sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
    
    outSE <- GLMPCA(se, L = 2, fam = "poi")
    
    expect_s3_class(metadata(outSE)$glmpca, "glmpca")
    gnames<-c("factors","loadings","coefX","coefZ","dev","glmpca_family") 
    expect_true(all(gnames %in% names(metadata(outSE)$glmpca)))
    
    outSCE <- GLMPCA(sce, L = 2, fam = "poi")
    
    expect_s3_class(metadata(outSE)$glmpca, "glmpca")
    expect_true("GLMPCA" %in% reducedDimNames(outSCE))
    
    outMAT <- GLMPCA(u, L = 2, fam = "poi")
    
    expect_s3_class(outMAT, "glmpca")
    expect_true(all(gnames %in% names(outMAT) ))
})
