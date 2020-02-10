#' Sum of vector elements.
#'
#' \code{sum} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.

suppressPackageStartupMessages(library(SingleCellExperiment))
library(Seurat)
source("./util/functions.R") #needed for compute_gene_info

#' Feature Selection by Deviance
#' 
#' In a typical single cell analysis, many of the features (genes) may have not be informative
#' about differences between cells. Feature selection seeks to identify which genes 
#' are the most informative. We define an informative feature as one that is poorly fit
#' by a multinomial model of constant expression across cells. We compute a deviance statistic
#' for each gene, then rank the genes in decreasing order. Genes with large deviances are
#' likely to be informative, and genes with small deviances can be discarded to speed up downstream
#' analyses such as dimension reduction.
#' 
#' @param object an object inheriting from \code{\link{SummarizedExperiment}} (such as \code{\link{SingleCellExperiment}}).
#' @param assay a character or integer specifying which assay to use for computing deviances of genes
#'   (default = 1).
#' @param fam a character specifying the model type to be used for calculating the
#'   deviance. Default is "binomial" which is the true marginal from a multinomial null model.
#'   The "poisson" alternative is an approximation to binomial.
#'   Option "geometric" is an overdispersed alternative that approximates a Dirichlet-multinomial null model.
#' @param nkeep how many informative genes should be retained? Default: all genes are retained.
#' @param ret should the result be returned as a SummarizedExperiment object or a data.frame?
#'   
#' @return If ret="sce", a subset of the original object containing only the informative genes, 
#'   and with deviance statistics appended to the rowData. 
#'   If ret="ranks", a data frame is returned with the deviance statistics.
#'   In both cases, the rows of the returned object will be resorted in decreasing deviance,
#'   so that the most informative genes are at the top.
#'    
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @export
featureSelectDeviance

compute_gene_info<-function(m,gmeta=NULL,fam=c("binomial","poisson","geometric")){
  #m a data matrix with genes=rows
  #gmeta a pre-existing data frame with gene-level metadata
  fam<-match.arg(fam)
  if(!is.null(gmeta)){ stopifnot(nrow(m)==nrow(gmeta)) }
  gnz<-Matrix::rowSums(m>0)
  sz<-compute_size_factors(m,fam)
  gof<-function(g){ gof_func(m[g,],sz,fam) }
  gof<-as.data.frame(t(vapply(1:nrow(m),gof,FUN.VALUE=rep(0.0,2))))
  #colnames(gof)<-c("deviance","pval")
  gmu<-Matrix::rowMeans(m)
  gvar<-apply(m,1,var)
  gfano<-ifelse(gvar>0 & gmu>0, gvar/gmu, 0)
  res<-cbind(nzsum=gnz,fano=gfano,gof)
  res$pval_fdr<-p.adjust(res$pval,"BH")
  if(is.null(gmeta)){ return(res) } else { return(cbind(gmeta,res)) }
}


filterDev<-function(sce,nkeep=nrow(sce),fam=c("binomial","poisson","geometric"),ret=c("sce","ranks")){
  fam<-match.arg(fam)
  ret<-match.arg(ret)
  gm<-compute_gene_info(counts(sce),gmeta=rowData(sce),fam=fam)
  o<-order(gm$deviance,decreasing=TRUE,na.last=FALSE)[1:nkeep]
  #NA deviance => badly fitting null model=> highly variable gene
  if(ret=="sce"){
    res<-sce[o,]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

