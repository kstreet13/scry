context("Test deviance feature selection")
require(SingleCellExperiment)
set.seed(1234)

test_that("featureSelection works", {
    ncells <- 100
    u <- matrix(rpois(2000, 3), ncol=ncells)
    de_genes<-5:10
    rownames(u)<-LETTERS[1:20]
    #create some differentially expressed genes
    u[de_genes,1:30]<-3*u[de_genes,1:30]
    u[de_genes,31:60]<-round(.5*u[de_genes,31:60])
    v <- log2(u + 1)
    m <- Matrix::Matrix(u,sparse=TRUE)
    se <- SummarizedExperiment(assays=list(logcounts=v,counts=u))
    sce <- SingleCellExperiment(assays=list(logcounts=v,counts=u,
                                            sparse_counts=m))
    h5<-as(m,"HDF5Matrix") #depends on package HDF5Matrix

    #check no error with different object types
    outU<-devianceFeatureSelection(u)
    expect_is(outU,"numeric")
    #all DE genes should have higher deviance than non-DE genes.
    expect_gt(min(outU[de_genes]), max(outU[-de_genes]))
    outM<-devianceFeatureSelection(m)
    outH5<-devianceFeatureSelection(h5)
    outSE<-devianceFeatureSelection(se,assay="counts")
    expect_is(outSE,"SummarizedExperiment")
    
    #check all object inputs give same output
    expect_equivalent(outM,outU)
    expect_equivalent(outH5,outU)
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
    #verify handling when nkeep > nrow(sce)
    nkeep <- nrow(sce)+10
    outSCE<-devianceFeatureSelection(sce,nkeep=nkeep)
    expect_equal(nrow(outSCE), nrow(sce))
    
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
    
    #check handling of batch argument
    expect_error(devianceFeatureSelection(sce, batch = rep(1, ncol(sce))),
                 ') is not TRUE')
    outSCE <- devianceFeatureSelection(sce, batch = factor(rep(1:2, length.out = ncol(sce))))
    expect_true("binomial_deviance" %in% names(rowData(outSCE)))
    
})