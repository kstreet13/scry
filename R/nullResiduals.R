# Helper functions for nullResiduals
.null_poisson_deviance_residuals_delayed <- function(m){
    # adapted from compute_size_factors
    lsz <- log(DelayedArray::colSums(m))
    sz <- exp(lsz-mean(lsz))
    
    # adapted from .null_residuals
    lambdahat <- DelayedArray::rowSums(m) / sum(sz)
    mhat <- DelayedArray(matrix(lambdahat)) %*%
        DelayedArray(matrix(sz, nrow = 1))
    
    # adapted from .poisson_deviance_residuals
    term1 <- m * log(m / mhat)
    term1[is.na(term1)] <- 0 #0*log(0)=0
    s2 <- 2 * (term1 - (m - mhat))
    sign(m - mhat) * sqrt(abs(s2))
}

.binomial_deviance_residuals <- function(X, p, n){
    #X a matrix, n is vector of length ncol(X)
    stopifnot(length(n) == ncol(X))
    if(is.matrix(X)){ #X is a dense matrix
        if(length(p) == nrow(X)){ 
            #p is a vector, length must match nrow(X)
            mu <- outer(p, n)
        } else if(!is.null(dim(p)) && dim(p) == dim(X)){
            # p is matrix, must have same dims as X
            mu <- t(t(p)*n)
        } else { 
            stop("dimensions of p and X must match!") 
        }
        term1 <- X*log(X/mu)
        term1[is.na(term1)] <- 0 #0*log(0)=0
        nx <- t(n - t(X))
        term2 <- nx*log(nx/outer(1-p, n))
        #this next line would only matter if all counts
        #were from a single gene, so not checking saves time.
        # term2[is.na(term2)] <- 0 
        return(sign(X-mu)*sqrt(2*(term1+term2)))
    } else { #X is a sparse Matrix or delayed Array
        stopifnot(length(p) == nrow(X))
        dr_func<-function(j){
            #j is a row index of X
            x <- X[j,]
            mu <- n*p[j]
            term1 <- x*log(x/mu)
            term1[is.na(term1)] <- 0 #0*log(0)=0
            nx <- n-x
            term2 <- nx*log(nx/(n*(1-p[j])))
            #this next line would only matter if all counts
            #were from a single gene, so not checking saves time.
            # term2[is.na(term2)] <- 0 
            sign(x-mu)*sqrt(2*(term1+term2))
        }
        #up to this point no large dense objects have been created
        #the last line here is to be modified to write each row
        #to a disk-based delayed Array object to avoid creating the big
        #dense object.
        return(t(vapply(seq_len(nrow(X)),dr_func,FUN.VALUE=0.0*n)))
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

.null_residuals <- function(m, fam = c("binomial", "poisson"), 
                           type = c("deviance", "pearson")){
    #m is a matrix, sparse Matrix, or delayed Array
    fam <- match.arg(fam); type <- match.arg(type)
    sz <- compute_size_factors(m, fam)
    if(fam == "binomial") {
        phat <- rowSums(m)/sum(sz)
        if(type == "deviance"){
            return(.binomial_deviance_residuals(m, phat, sz))
        } else { #pearson residuals
            if(is.matrix(m)){
                mhat <- outer(phat, sz)
                res <- (m-mhat)/sqrt(mhat*(1-phat))
                res[is.na(res)] <- 0 #case of 0/0
                return(res)
            } else { #if m is sparse Matrix or delayed Array
                pr_func<-function(j){
                    mhat <- phat[j]*sz
                    res <- (m[j,]-mhat)/sqrt(mhat*(1-phat[j]))
                    res[is.na(res)] <- 0 #case of 0/0
                }
                #this last line is the only part where the full dense object
                #is instantiated in memory, replace with row-by-row write to 
                #delayedArray
                return(t(vapply(seq_len(nrow(m),pr_func,FUN.VALUE=0.0*sz))))
            } #end sparse binomial Pearson block
        } #end general binomial Pearson residuals block
    } else { #fam == "poisson"
        lambda<-rowSums(m)/sum(sz)
        if(is.matrix(m)){ #dense data matrix
            mhat <- outer(lambda, sz) 
            if(type == "deviance"){ 
                return(.poisson_deviance_residuals(m, mhat))
            } else { #pearson residuals
                res <- (m-mhat)/sqrt(mhat)
                res[is.na(res)] <- 0 #case of 0/0
                return(res)
            } #end dense Poisson Pearson residuals block
        } else { #case where m is a sparse Matrix or delayed Array
            lambda<- rowSums(m)/sum(sz)
            if(type == "deviance"){
                rfunc<-function(j){
                    .poisson_deviance_residuals(m[j,], lambda[j]*sz)
                }
            } else { #pearson residuals
                rfunc<-function(j){
                    mhat<- lambda[j]*sz
                    res<- (m[j,]-mhat)/sqrt(mhat)
                    res[is.na(res)] <- 0
                    res
                }
            }
            #up to this point no dense objects created in memory, 
            #modify below line
            #to write each row to a disk based delayedArray
            return(t(vapply(seq_len(nrow(m)),rfunc,FUN.VALUE=0.0*sz)))
        } #end sparse Poisson block
    } #end general Poisson block
}

.null_residuals_batch <- function(m, fam=c("binomial", "poisson"),
                                 type=c("deviance", "pearson"), batch=NULL){
    #null residuals but with batch indicator (batch=a factor)
    fam <- match.arg(fam); type <- match.arg(type)
    if(is.null(batch)){
        return(.null_residuals(m, fam = fam, type = type))
    } else { #case where there is more than one batch
        stopifnot(length(batch) == ncol(m) && is(batch, "factor"))
        res <- matrix(0.0, nrow = nrow(m), ncol = ncol(m))
        for(b in levels(batch)){
            idx <- (batch == b)
            res[, idx] <- .null_residuals(m[, idx], fam = fam, type = type)
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
#'   genes int he rows and samples in the columns. Specialized methods are 
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
          definition = function(object, assay = 1,
                                fam = c("binomial", "poisson"), 
                                type = c("deviance", "pearson"), 
                                batch = NULL){
              fam <- match.arg(fam); type <- match.arg(type)
              # m <- as.matrix(assay(object, assay))
              name <- paste(fam, type, "residuals", sep="_")
              assay(object, name) <- .null_residuals_batch(m, fam, type, batch)
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
              #.null_residuals_batch(as.matrix(object), fam, type, batch)
              .null_residuals_batch(object, fam, type, batch)
          })

#' @rdname nullResiduals
#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setMethod(f = "nullResiduals", 
          signature = signature(object = "ANY"), 
          definition = function(object, fam = c("binomial", "poisson"), 
                                type = c("deviance", "pearson"), 
                                batch = NULL){
              if(!is(x, "HDF5Matrix") && !is(x, "DelayedMatrix")) {
                  stop("x is of type ", class(x), ", currently not supported")
                } else {
                    fam <- match.arg(fam); type <- match.arg(type)
                    .null_residuals_batch(object, fam, type, batch)
                }
          })