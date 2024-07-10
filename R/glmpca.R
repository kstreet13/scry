#' @title Generalized principal component analysis
#' @rdname GLMPCA
#' @description This function implements the GLM-PCA dimensionality reduction
#'   method for high-dimensional count data. This is a wrapper for
#'   \code{\link[glmpca]{glmpca}}.
#' 
#' @param object A \code{\link{SingleCellExperiment}} or
#'   \link{SummarizedExperiment} object. Alternatively, a matrix-like object
#'   of non-negative integer counts (such as a sparse \code{\link{Matrix}}).
#' @param L the desired number of latent dimensions (integer).
#' @param assay a character or integer specifying which assay to use for GLM-PCA
#'   (default = 'counts'). Ignored if \code{object} is a matrix.
#' @param ... further arguments passed to \code{\link[glmpca]{glmpca}}
#'
#' @return The original \code{SingleCellExperiment} or
#'   \code{SummarizedExperiment} object with the GLM-PCA results added to the
#'   \code{metadata} slot. If the original input was a
#'   \code{SingleCellExperiment}, then a new \code{reducedDim} element called
#'   \code{"GLMPCA"} will be added, representing the GLM-PCA \code{factors}. If
#'   the input was a matrix, output matches that of
#'   \code{\link[glmpca]{glmpca}}.
#'   
#' @examples 
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=u))
#' GLMPCA(sce, L = 2)
#' 
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom methods is
#' @export
setMethod(f = "GLMPCA","SummarizedExperiment",
          definition = function(object, L, assay = "counts", ...){
              res <- glmpca::glmpca(Y = assay(object, assay), L = L, ...)
              if(is(object, 'SingleCellExperiment')){
                  reducedDim(object, "GLMPCA") <- as.matrix(res$factors)
              }
              object@metadata$glmpca <- res
              return(object)
          })

#this function is only needed because the glmpca signature has "Y" not "object"
glmpca_wrapper<-function(object, L, ...){
    glmpca::glmpca(object, L, ...)
}

#' @rdname GLMPCA
#' @export
setMethod(f = "GLMPCA", "matrix", definition=glmpca_wrapper)

#' @rdname GLMPCA
#' @export
setMethod(f = "GLMPCA", "Matrix", definition=glmpca_wrapper)
