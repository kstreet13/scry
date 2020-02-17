context("Test devianceFeatureSelection, GLMPCA, and nullResiduals")
set.seed(1234)

test_that("featureSelection works with intended input types", {
  ncells <- 100
  u <- matrix(rpois(2000, 3), ncol=ncells)
  de_genes<-5:10
  rownames(u)<-LETTERS[1:20]
  #create some differentially expressed genes
  u[de_genes,1:30]<-3*u[de_genes,1:30]
  u[de_genes,31:60]<-round(.5*u[de_genes,31:60])
  v <- log2(u + 1)
  require(Matrix)
  m <- Matrix(u)
  require(SingleCellExperiment)
  se <- SummarizedExperiment(assays=list(logcounts=v,counts=u))
  sce <- SingleCellExperiment(assays=list(logcounts=v,counts=u,sparse_counts=m))
  
  #check no error with different object types
  outU<-devianceFeatureSelection(u)
  expect_is(outU,"numeric")
  #all DE genes should have higher deviance than non-DE genes.
  expect_gt(min(outU[de_genes]), max(outU[-de_genes]))
  outM<-devianceFeatureSelection(m)
  
  outSE<-devianceFeatureSelection(se,assay="counts")
  expect_is(outSE,"SummarizedExperiment")
  
  #check all object inputs give same output
  expect_equivalent(outM,outU)
  expect_equivalent(rowData(outSE)[,"binomial_deviance"], outU)
  
  #check different input parameters don't throw errors
  outSCE1 <- devianceFeatureSelection(sce, assay="counts", fam = "poisson")
  outSCE2 <- devianceFeatureSelection(sce, assay=3, fam = "poisson") #sparse counts
  #check sparse counts and counts give same result
  expect_equivalent(rowData(outSCE1)[,"poisson_deviance"],rowData(outSCE2)[,"poisson_deviance"])
  outSCE <- devianceFeatureSelection(sce, assay="sparse_counts", fam = "binomial")
  
  #check SCE output equivalent to matrix output
  expect_equivalent(rowData(outSCE)[,"binomial_deviance"],outU)
  
  #check it doesn't alter the properties of the original SCE object.
  expect_equal(rownames(outSCE),rownames(sce))
  expect_equal(assayNames(outSCE),assayNames(sce))
  
  #check sorting and subsetting works properly
  #sorting but no subsetting
  outSCE<-devianceFeatureSelection(sce,assay="sparse_counts",sorted=TRUE)
  #verify de_genes have top deviance
  expect_true(all(rownames(u)[de_genes] %in% rownames(outSCE)[1:length(de_genes)]))
  #verify it's in decreasing order of deviance
  outSCE_devs<-rowData(outSCE)[,"binomial_deviance"]
  expect_equivalent(outSCE_devs,sort(outSCE_devs,decreasing=TRUE)) 
  
  #subsetting leading to automatic sorting
  nkeep<-length(de_genes)+2
  outSCE<-devianceFeatureSelection(sce,assay="sparse_counts",nkeep=nkeep)
  expect_equal(nrow(outSCE),nkeep)
  #verify de_genes have top deviance
  expect_true(all(rownames(u)[de_genes] %in% rownames(outSCE)[1:length(de_genes)]))
  #verify it's in decreasing order of deviance
  outSCE_devs<-rowData(outSCE)[,"binomial_deviance"]
  expect_equivalent(outSCE_devs,sort(outSCE_devs,decreasing=TRUE))
  
  #introduce batch effect- this is not working yet...
  # batch_genes<-1:3 #make sure this doesn't overlap with de_genes
  # batch<-rep("b",ncells)
  # odds<-seq_along(batch)%%2==1
  # batch[odds]<-"a"
  # u[batch_genes,odds]<-2*u[batch_genes,odds] #batch effect
  # sce <- SingleCellExperiment(assays=list(counts=Matrix(u)))
  # outSCE<-devianceFeatureSelection(sce,nkeep=length(de_genes))
  # #verify batch effect genes appear as high deviance without adjustment
  # expect_true(any(rownames(u)[batch_genes] %in% rownames(outSCE)))
  # outSCE<-devianceFeatureSelection(sce,nkeep=length(de_genes),batch=batch)
  # #verify batch effect genes do not appear as high deviance with adjustment
  # expect_false(any(rownames(u)[batch_genes] %in% rownames(outSCE)))
})

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

test_that("Cell type prediction works with intended input types", {
	# make/import toy data
	
	# run predictType on toy data
	
	# check validity of outputs
	
})
