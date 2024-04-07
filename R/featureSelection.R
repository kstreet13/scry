#compute deviance for each gene (row) in a matrix-like object...
#...or SummarizedExperiment

#' @importFrom methods as
#' @importFrom Matrix Diagonal
sparseBinomialDeviance<-function(X,sz){
    #X has features in cols, observations in rows
    #assume X is a sparseMatrix object
    X<-as(X,"CsparseMatrix")
    LP<-L1P<-Matrix::Diagonal(x = 1/sz) %*% X #recycling
    LP@x<-log(LP@x) #log transform nonzero elements only
    L1P@x<-log1p(-L1P@x) #rare case: -Inf if only a single gene nonzero in a cell
    ll_sat<-Matrix::colSums(X*(LP-L1P)+sz*L1P, na.rm=TRUE)
    sz_sum<-sum(sz)
    feature_sums<-Matrix::colSums(X)
    p<-feature_sums/sz_sum
    l1p<-log1p(-p)
    ll_null<-feature_sums*(log(p)-l1p)+sz_sum*l1p
    2*(ll_sat-ll_null)
}

denseBinomialDeviance<-function(X,sz){
    #X has features in cols, observations in rows
    P<-X/sz
    L1P<-log1p(-P)
    ll_sat<-DelayedArray::colSums(X*(log(P)-L1P)+sz*L1P, na.rm=TRUE)
    sz_sum<-sum(sz)
    feature_sums<-DelayedArray::colSums(X)
    p<-feature_sums/sz_sum
    l1p<-log1p(-p)
    ll_null<-feature_sums*(log(p)-l1p)+sz_sum*l1p
    2*(ll_sat-ll_null)
}

#' @importFrom methods as
sparsePoissonDeviance<-function(X,sz){
    #X has features in cols, observations in rows
    X<-as(X,"CsparseMatrix")
    LP<-X/sz #recycling
    LP@x<-log(LP@x) #log transform nonzero elements only
    ll_sat<-Matrix::colSums(X*LP, na.rm=TRUE)
    feature_sums<-Matrix::colSums(X)
    ll_null<-feature_sums*log(feature_sums/sum(sz))
    2*(ll_sat-ll_null)
}

densePoissonDeviance<-function(X,sz){
    #X has features in cols, observations in rows
    ll_sat<-DelayedArray::colSums(X*log(X/sz), na.rm=TRUE)
    feature_sums<-DelayedArray::colSums(X)
    ll_null<-feature_sums*log(feature_sums/sum(sz))
    2*(ll_sat-ll_null)
}

#' @importFrom Matrix t
#' @importFrom methods is
.compute_deviance<-function(m,fam=c("binomial","poisson")){
    #m is either a Matrix or matrix object (later: support DelayedArrays)
    #m a data matrix with genes=rows
    fam <- match.arg(fam)
    sz <- colSums(m)
    if(fam=="poisson"){
        lsz<-log(sz)
        #make geometric mean of sz be 1 for poisson
        sz <- exp(lsz-mean(lsz))
    }
    m<-t(m) #column slicing faster than row slicing for matrix in R.
    #note: genes are now in the COLUMNS of m
    if(is(m,"sparseMatrix")){
        if(fam=="binomial"){
            out <- sparseBinomialDeviance(m,sz)
        } else { #fam=="poisson"
            out <- sparsePoissonDeviance(m,sz)
        }
    } else { #m is either 1) an ordinary dense array or matrix
        # 2) a non-sparse Matrix object like dgeMatrix
        # 3) a dense object like HDF5Array (on disk) or DelayedArray (in memory)
        if(fam=="binomial"){
            out <- denseBinomialDeviance(m,sz)
        } else { #fam=="poisson"
            out <- densePoissonDeviance(m,sz)
        }
    }
    out[which(is.na(out))] <- 0
    return(out)
}

.compute_deviance_batch<-function(object,
                                  fam=c("binomial","poisson"),
                                  batch=NULL){
    #deviance but with batch indicator (batch=a factor)
    m<-object; rm(object)
    fam<-match.arg(fam)
    stopifnot(is.null(batch) || is(batch,"factor"))
    if(is.null(batch)){
        return(.compute_deviance(m,fam=fam))
    } else { #case where there is more than one batch
        stopifnot(length(batch)==ncol(m))
        #each row is a gene, each column is deviance within a batch.
        res<-matrix(0.0,nrow=nrow(m),ncol=nlevels(batch))
        for(i in seq_along(levels(batch))){
            b<-levels(batch)[i]
            res[,i]<-.compute_deviance(m[,batch==b],fam=fam)
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
#' @param object an object inheriting from \code{\link{SummarizedExperiment}}
#'   (such as
#'   \code{\link{SingleCellExperiment}}). Alternatively, a matrix or matrix-like
#'   object (such as a sparse \code{\link{Matrix}}) of non-negative integer
#'   counts.
#' @param assay a string or integer specifying which assay contains the count
#'   data (default = 'counts'). Ignored if \code{object} is a matrix-like
#'   object.
#' @param fam a string specifying the model type to be used for calculating the
#'   residuals. Binomial (the default) is the closest approximation to
#'   multinomial, but Poisson may be faster to compute and often is very similar
#'   to binomial.
#' @param batch an optional factor indicating batch membership of observations.
#'   If provided, the null model is computed within each batch separately to
#'   regress out the batch effect from the resulting deviance statistics.
#' @param nkeep integer, how many informative features should be retained?
#'   Default: all features are retained if set to NULL. Ignored if \code{object}
#'   is a matrix-like object.
#' @param sorted logical, should the \code{object} be returned with rows sorted
#'   in decreasing order of deviance? Default: FALSE, unless nkeep is specified,
#'   in which case it is forced to be TRUE. Ignored for matrix-like inputs.
#'
#' @return The original \code{SingleCellExperiment} or
#'   \code{SummarizedExperiment} object with the deviance statistics for each
#'   feature appended to the rowData. The new column name will be either
#'   binomial_deviance or poisson_deviance. If the input was a matrix-like
#'   object, output is a numeric vector containing the deviance statistics for
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
setMethod("devianceFeatureSelection", "SummarizedExperiment",
          definition = function(object, assay = "counts",
                                fam = c("binomial", "poisson"), batch = NULL,
                                nkeep = NULL, sorted = FALSE){
              fam<-match.arg(fam)
              m <- assay(object, assay)
              dev<-.compute_deviance_batch(m, fam, batch)
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
setMethod("devianceFeatureSelection", "matrix",
          definition=.compute_deviance_batch)

#' @rdname devianceFeatureSelection
#' @export
setMethod("devianceFeatureSelection", "Matrix",
          definition=.compute_deviance_batch)

#' @rdname devianceFeatureSelection
#' @export
setMethod("devianceFeatureSelection", "DelayedArray",
          definition=.compute_deviance_batch)
