poisson_deviance<-function(x,mu,sz){
    #assumes log link and size factor sz on the same scale as x (not logged)
    #stopifnot(all(x>=0 & sz>0))
    2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

binomial_deviance<-function(x,p,n){
    term1<-sum(x*log(x/(n*p)), na.rm=TRUE)
    nx<-n-x
    term2<-sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
    2*(term1+term2)
}

#' @importFrom Matrix t
compute_deviance<-function(m,fam=c("binomial","poisson")){
    #m a data matrix with genes=rows
    fam<-match.arg(fam)
    sz<-compute_size_factors(m,fam)
    sz_sum<-sum(sz)
    m<-Matrix::t(m) #column slicing faster than row slicing for matrix in R.
    #note: genes are now in the COLUMNS of m
    dev_binom<-function(g){
        x<-m[,g]
        binomial_deviance(x,sum(x)/sz_sum,sz)
    }
    dev_poi<-function(g){
        x<-m[,g]
        poisson_deviance(x,sum(x)/sz_sum,sz)
    }
    dev<-if(fam=="binomial"){ dev_binom } else { dev_poi }
    #we can't just use apply() b/c it will coerce a sparse Matrix to dense.
    #note: would be possible to parallelize over cols of m here
    #computation for each gene is independent. Future enhancement?
    vapply(seq_len(ncol(m)),dev,FUN.VALUE=0.0) #numeric vector
}

compute_deviance_batch<-function(m,fam=c("binomial","poisson"),batch=NULL){
    #deviance but with batch indicator (batch=a factor)
    fam<-match.arg(fam)
    stopifnot(is.null(batch) || is(batch,"factor"))
    if(is.null(batch)){
        return(compute_deviance(m,fam=fam))
    } else { #case where there is more than one batch
        stopifnot(length(batch)==ncol(m))
        #each row is a gene, each column is deviance within a batch.
        res<-matrix(0.0,nrow=nrow(m),ncol=nlevels(batch))
        for(i in seq_along(levels(batch))){
            b<-levels(batch)[i]
            res[,i]<-compute_deviance(m[,batch==b],fam=fam)
        }
        #deviance is additive across subsets of observations
        return(rowSums(res)) 
    }
}

#' @title Feature selection by approximate multinomial deviance
#' @rdname devianceFeatureSelection
#' @description Computes a deviance statistic for each row feature (such as a
#'   gene) for count data based on a multinomial null model that assumes each
#'   feature has a constant rate. Features with large deviance are likely to be
#'   informative. Uninformative, low deviance features can be discarded to speed
#'   up downstream analyses and reduce memory footprint.
#' 
#' @param object an object inheriting from \link{SummarizedExperiment} (such as
#'   \code{\link{SingleCellExperiment}}). Alternatively, a matrix or sparse
#'   Matrix of integer counts.
#' @param assay a string or integer specifying which assay contains the count
#'   data (default = 1). Ignored if \code{object} is a matrix or sparse Matrix.
#' @param fam a string specifying the model type to be used for calculating the
#'   residuals. Binomial (the default) is the closest approximation to
#'   multinomial, but Poisson may be faster to compute and often is very similar
#'   to binomial.
#' @param batch an optional factor indicating batch membership of observations.
#'   If provided, the null model is computed within each batch separately to
#'   regress out the batch effect from the resulting deviance statistics.
#' @param nkeep integer, how many informative features should be retained?
#'   Default: all features are retained if set to NULL. Ignored if \code{object}
#'   is a matrix or sparse Matrix.
#' @param sorted logical, should the \code{object} be returned with rows sorted
#'   in decreasing order of deviance? Default: FALSE, unless nkeep is specified,
#'   in which case it is forced to be TRUE. Ignored for matrix and sparse Matrix
#'   inputs.
#'   
#' @return The original \code{SingleCellExperiment} or
#'   \code{SummarizedExperiment} object with the deviance statistics for each
#'   feature appended to the rowData. The new column name will be either
#'   binomial_deviance or poisson_deviance. If the input was a matrix or sparse
#'   Matrix, output is a numeric vector containing the deviance statistics for
#'   each row.
#'
#' @details In a typical single-cell analysis, many of the features (genes) may
#'   not be informative about differences between observations (cells). Feature
#'   selection seeks to identify which genes are the most informative. We define
#'   an informative gene as one that is poorly fit by a multinomial model of
#'   constant expression across cells within each batch. We compute a deviance
#'   statistic for each gene. Genes with high deviance are more informative.
#'   
#' @references 
#' Townes FW, Hicks SC, Aryee MJ, and Irizarry RA (2019). Feature
#' Selection and Dimension Reduction for Single Cell RNA-Seq based on a
#' Multinomial Model. \emph{Genome Biology}
#' \url{https://doi.org/10.1186/s13059-019-1861-6}
#' 
#' @examples 
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=u))
#' devianceFeatureSelection(sce)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @export
setMethod(f = "devianceFeatureSelection",
          signature = signature(object = "SummarizedExperiment"),
          definition = function(object, assay = 1, 
                                fam = c("binomial", "poisson"), batch = NULL, 
                                nkeep = NULL, sorted = FALSE){
              fam<-match.arg(fam)
              m <- assay(object, assay)
              dev<-compute_deviance_batch(m, fam, batch)
              name<-paste(fam, "deviance", sep="_")
              rd<-rowData(object)
              rd[name]<-dev
              rowData(object)<-rd
              if(!is.null(nkeep) && nkeep>=length(dev)){
                  nkeep<-NULL
              } #user wants to keep all features
              if(!is.null(nkeep)){ sorted<-TRUE } #force sorting if we are 
              # taking a subset of rows
              if(sorted){ 
                  o<-order(dev,decreasing=TRUE)
                  object<-object[o,]
                  if(is.null(nkeep)){ #sorting but no subsetting
                      return(object) 
                  } else { #sorting and subsetting
                      return(object[seq_len(nkeep),]) 
                  }
              } else { #no sorting, no subsetting
                  return(object)
              }
          })

#' @rdname devianceFeatureSelection
#' @export
setMethod(f = "devianceFeatureSelection",
          signature = signature(object = "matrix"),
          definition = function(object, fam = c("binomial", "poisson"),
                                batch = NULL){
              fam<-match.arg(fam)
              compute_deviance_batch(object,fam,batch)
          })

#' @rdname devianceFeatureSelection
#' @export
setMethod(f = "devianceFeatureSelection",
          signature = signature(object = "Matrix"),
          definition = function(object, fam = c("binomial", "poisson"),
                                batch = NULL){
              fam<-match.arg(fam)
              compute_deviance_batch(object,fam,batch)
          })