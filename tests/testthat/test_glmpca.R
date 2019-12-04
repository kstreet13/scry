context("Test GLMPCA and devianceResiduals")
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



test_that("devianceResiduals works with intended input types", {
	ncells <- 100
	u <- matrix(rpois(2000, 5), ncol=ncells)
	v <- log2(u + 1)
	
	require(SummarizedExperiment)
	se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
	sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
	
	outSE <- devianceResiduals(se, fam = "bin")
	expect_true(all(c('deviance','dev_pval','dev_qval') %in% 
						names(colData(outSE))))

	expect_message(outSE <- devianceResiduals(se, assay = 'logcounts', fam = "poi",
							   recalcSizeFactors = FALSE), 
				   'Cannot use existing size factors')
	expect_true(all(c('deviance','dev_pval','dev_qval') %in% 
						names(colData(outSE))))
	
	outSE <- devianceResiduals(se, assay = 2, fam = "geo")
	expect_true(all(c('deviance','dev_pval','dev_qval') %in% 
						names(colData(outSE))))
	

	# pathological cases
	
})