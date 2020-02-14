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
})

