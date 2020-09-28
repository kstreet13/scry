# Helper functions for nullResiduals
.binomial_deviance_residuals <- function(X, p, n){
    #X a matrix, n is vector of length ncol(X)
    stopifnot(length(n) == ncol(X))
    stopifnot(length(p) == nrow(X))
    #p is a vector, length must match nrow(X)
    if(is.matrix(X) | is(X,"Matrix")){ #X is in-memory
        mu <- outer(p, n)
        term1 <- X*log(X/mu)
        term1[is.na(term1)] <- 0 #0*log(0)=0
        nx <- t(n - t(X))
        term2 <- nx*log(nx/outer(1-p, n))
        #this next line would only matter if all counts
        #were from a single gene, so not checking saves time.
        # term2[is.na(term2)] <- 0
        res <- sign(X-mu)*sqrt(2*(term1+term2))
        res[is.na(res)] <- 0 # handle cases where term1+term2 = -epsilon
        return(res)
    } else { #X delayed Array / out-of-memory
        mu <- BiocSingular::LowRankMatrix(
            DelayedArray(matrix(p)),
            DelayedArray(matrix(n)))
        term1 <- X*log(X/mu)
        term1[is.na(term1)] <- 0 #0*log(0)=0
        nx <- t(n - t(X))
        term2 <- nx*log(nx/BiocSingular::LowRankMatrix(
            DelayedArray(matrix(1-p)),
            DelayedArray(matrix(n))))
        #this next line would only matter if all counts
        #were from a single gene, so not checking saves time.
        term2[is.na(term2)] <- 0 
        res <- sign(X-mu)*sqrt(2*(term1+term2))
        res[is.na(res)] <- 0 # handle cases where term1+term2 = -epsilon
        return(res)
    }
}

.poisson_deviance_residuals <- function(x, xhat){
    #x, xhat assumed to be same dimension
    #sz <- exp(offsets)
    #xhat <- lambda*sz
    term1 <- x*log(x/xhat)
    term1[is.na(term1)] <- 0 #0*log(0)=0
    s2 <- 2*(term1-(x-xhat))
    sign(x-xhat)*sqrt(abs(s2))
}

#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom DelayedArray colSums
.null_residuals <- function(m, fam = c("binomial", "poisson"),
                            type = c("deviance", "pearson"),
                            size_factors=NULL){
    #m is a matrix or sparse Matrix
    fam <- match.arg(fam); type <- match.arg(type)
    if(is.null(size_factors)) {
        sz <- colSums(m)
    } else {
        sz <- size_factors
    }

    if(fam == "binomial") {
        phat <- rowSums(m)/sum(sz)
        if(type == "deviance"){
            return(.binomial_deviance_residuals(m, phat, sz))
        } else { #pearson residuals
            # make mhat
            if(is.matrix(m) | is(m,"Matrix")){
                mhat <- outer(phat, sz)
            } else { #if m is delayed Array / out-of-memory
                mhat <- BiocSingular::LowRankMatrix(
                    DelayedArray(matrix(phat)),
                    DelayedArray(matrix(sz)))
            }
            res <- (m-mhat)/sqrt(mhat*(1-phat))
            res[is.na(res)] <- 0 #case of 0/0
            return(res)
        }
    } else { #fam == "poisson"
        lsz <- log(sz)
        #make geometric mean of sz be 1 for poisson
        sz <- exp(lsz-mean(lsz))
        lambda <- rowSums(m) / sum(sz)
        # make mhat
        if(is.matrix(m) | is(m,"Matrix")){ #dense data matrix
            mhat <- outer(lambda, sz)
        } else { #case where m is delayed Array
            mhat <- BiocSingular::LowRankMatrix(
                DelayedArray(matrix(lambda)),
                DelayedArray(matrix(sz)))
        }
        if(type == "deviance"){
            return(.poisson_deviance_residuals(m, mhat))
        } else { # pearson residuals
            res <- (m-mhat)/sqrt(mhat)
            res[is.na(res)] <- 0 #case of 0/0
            return(res)
        }
    } #end general Poisson block
}


.null_residuals_batch <- function(m, fam=c("binomial", "poisson"),
                                  type=c("deviance", "pearson"), batch=NULL, size_factors=NULL){
    #null residuals but with batch indicator (batch=a factor)
    fam <- match.arg(fam); type <- match.arg(type)
    if(is.null(batch)){
        return(.null_residuals(m, fam = fam, type = type, size_factors = size_factors))
    } else { #case where there is more than one batch
        stopifnot(length(batch) == ncol(m) && is(batch, "factor"))
        res <- matrix(0.0, nrow = nrow(m), ncol = ncol(m))
        for(b in levels(batch)){
            idx <- which(batch == b)
            res[, idx] <- .null_residuals(m[, idx], fam = fam, type = type, size_factors = size_factors[idx])
        }
        return(res)
    }
}

#' @title Residuals from an approximate multinomial null model
#' @rdname nullResiduals
#' @description Computes deviance or Pearson residuals for count data based on a
#'   multinomial null model that assumes each feature has a constant rate. The
#'   residuals matrix can be analyzed with standard PCA as a fast approximation
#'   to GLM-PCA.
#'
#' @param object The object on which to compute residuals. It can be a
#'   matrix-like object (e.g. matrix, Matrix, DelayedMatrix, HDF5Matrix) with
#'   genes in the rows and samples in the columns. Specialized methods are
#'   defined for objects inheriting from \link{SummarizedExperiment} (such as
#'   \code{\link{SingleCellExperiment}}).
#' @param assay a string or integer specifying which assay contains the count
#'   data (default = 1). Ignored if \code{object} is a matrix.
#' @param fam a string specifying the model type to be used for calculating the
#'   residuals. Binomial (the default) is the closest approximation to
#'   multinomial, but Poisson may be faster to compute and often is very similar
#'   to binomial.
#' @param type should deviance or Pearson residuals be used?
#' @param batch an optional factor indicating batch membership of observations.
#'   If provided, the null model is computed within each batch separately to
#'   regress out the batch effect from the resulting residuals.
#'
#' @return The original \code{SingleCellExperiment} or
#'   \code{SummarizedExperiment} object with the residuals appended as a new
#'   assay. The assay name will be fam_type_residuals (eg,
#'   binomial_deviance_residuals). If the input was a matrix, output is a dense
#'   matrix containing the residuals.
#'
#' @details This function should be used only on the un-normalized counts.
#'  It was originally designed for single-cell RNA-seq counts
#'  obtained by the use of unique molecular identifiers (UMIs) and has not been
#'  tested on read count data without UMIs or other data types.
#'
#'  Note that even though sparse Matrix objects are accepted as input,
#'  they are internally coerced to dense matrix before processing,
#'  because the output
#'  is always a dense matrix since the residuals transformation
#'  is not sparsity preserving.
#'  To avoid memory issues, it is recommended to perform feature selection first
#'  and subset the number of features to a smaller size prior to computing the
#'  residuals.
#'
#' @references Townes FW, Hicks SC, Aryee MJ, and Irizarry RA (2019). Feature
#' Selection and Dimension Reduction for Single Cell RNA-Seq based on a
#' Multinomial Model. \emph{Genome Biology}
#' \url{https://doi.org/10.1186/s13059-019-1861-6}
#'
#' @examples
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=u))
#' nullResiduals(sce)
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @export
setMethod(f = "nullResiduals",
          signature = signature(object = "SummarizedExperiment"),
          definition = function(object, assay = "counts",
                                fam = c("binomial", "poisson"),
                                type = c("deviance", "pearson"),
                                batch = NULL){
              fam <- match.arg(fam); type <- match.arg(type)
              name <- paste(fam, type, "residuals", sep="_")
              assay(object, name) <- .null_residuals_batch(assay(object, assay), fam, type, batch)
              object
          })

#' @rdname nullResiduals
#' @import SingleCellExperiment
#' @export
setMethod(f = "nullResiduals",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object, assay = "counts",
                                fam = c("binomial", "poisson"),
                                type = c("deviance", "pearson"),
                                batch = NULL){
              fam <- match.arg(fam); type <- match.arg(type)
              name <- paste(fam, type, "residuals", sep="_")
              tmp <- .null_residuals_batch(assay(object, assay), fam, type, batch, sizeFactors(object))
              rownames(tmp) <- rownames(object)
              colnames(tmp) <- colnames(object)
              assay(object, name) <- tmp
              object
          })

#' @rdname nullResiduals
#' @export
setMethod(f = "nullResiduals",
          signature = signature(object = "matrix"),
          definition = function(object, fam = c("binomial", "poisson"),
                                type = c("deviance", "pearson"),
                                batch = NULL){
              fam <- match.arg(fam); type <- match.arg(type)
              .null_residuals_batch(object, fam, type, batch)
          })

#' @rdname nullResiduals
#' @export
setMethod(f = "nullResiduals",
          signature = signature(object = "Matrix"),
          definition = function(object, fam = c("binomial", "poisson"),
                                type = c("deviance", "pearson"),
                                batch = NULL){
              fam <- match.arg(fam); type <- match.arg(type)

              warning("The matrix with the residuals will not be sparse.\n",
                      "Note that this will result in an object substantially larger than the input.\n",
                      "In some extreme cases this could result in an out-of-memory error.\n",
                      "We recommend using HDF5Matrix for large datasets. See the vignette for an example.")
              .null_residuals_batch(object, fam, type, batch)
          })

#' @rdname nullResiduals
#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importFrom DelayedArray DelayedArray
#' @importClassesFrom BiocSingular LowRankMatrix
setMethod(f = "nullResiduals",
          signature = signature(object = "ANY"),
          definition = function(object, fam = c("binomial", "poisson"),
                                type = c("deviance", "pearson"),
                                batch = NULL){
              if(!is(object, "HDF5Matrix") && !is(object, "DelayedMatrix")) {
                  stop("object is of type ", class(object),
                       ", currently not supported")
              } else {
                  fam <- match.arg(fam); type <- match.arg(type)
                  .null_residuals_batch(object, fam, type, batch)
              }
          })
