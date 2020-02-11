context("Test GLMPCA and nullResiduals")
set.seed(1234)

test_that("GLMPCA works with intended input types", {
	ncells <- 100
	u <- matrix(rpois(2000, 5), ncol=ncells)
	v <- log2(u + 1)
	
	require(SingleCellExperiment)
	se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
	sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
	
	outSE <- GLMPCA(se, L = 2, fam = "poi")
	
	expect_is(metadata(outSE)$glmpca, "list")
	expect_equal(names(metadata(outSE)$glmpca),
				 c("factors","loadings","coefX","coefZ","dev","family"))
	
	outSCE <- GLMPCA(sce, L = 2, fam = "poi")
	
	expect_is(metadata(outSE)$glmpca, "list")
	expect_true("GLMPCA" %in% reducedDimNames(outSCE))
	
	outMAT <- GLMPCA(u, L = 2, fam = "poi")
	
	expect_is(outMAT, "list")
	expect_equal(names(outMAT),
				 c("factors","loadings","coefX","coefZ","dev","family"))
	
	
	# pathological cases

})

test_that("nullResiduals works with intended input types", {
	ncells <- 100
	u <- matrix(rpois(2000, 5), ncol=ncells)
	v <- log2(u + 1)
	
	require(SummarizedExperiment)
	se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
	sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
	
	outSE <- nullResiduals(se, fam = "binomial")
	expect_true("binomial_deviance_residuals" %in% assayNames(outSE))
	outSE <- nullResiduals(se, fam = "binomial", type = "pearson")
	expect_true("binomial_pearson_residuals" %in% assayNames(outSE))
	outSE <- nullResiduals(se, fam = "poisson")
	expect_true("poisson_deviance_residuals" %in% assayNames(outSE))
	outSE <- nullResiduals(se, fam = "poisson", type = "pearson")
	expect_true("poisson_pearson_residuals" %in% assayNames(outSE))

})
