context("Test GLMPCA")
set.seed(1234)

test_that("glmpca works with intended input types", {
	ncells <- 100
	u <- matrix(rpois(2000, 5), ncol=ncells)
	v <- log2(u + 1)
	
	se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
	sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
	
	outSE <- GLMPCA(se, L = 2, fam = "poi")
	
	expect_is(metadata(outSE)$glmpca, "list")
	expect_equal(names(metadata(outSE)$glmpca),
				 c("factors","loadings","coefX","coefZ","dev","family"))
	
	outSCE <- GLMPCA(sce, L = 2, fam = "poi")
	
	expect_is(metadata(outSE)$glmpca, "list")
	expect_true("GLMPCA" %in% reducedDimNames(outSCE))
	
	# pathological cases

})