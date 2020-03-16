#' @title Generalized principal component analysis
#' @rdname GLMPCA
#' @description This function implements the GLM-PCA dimensionality reduction
#'   method for high-dimensional count data. This is a wrapper for
#'   \code{\link[glmpca]{glmpca}}.
#' 
#' @param object A \code{\link{SingleCellExperiment}} or
#'   \link{SummarizedExperiment} object. Alternatively, a matrix of integer
#'   counts.
#' @param assay a character or integer specifying which assay to use for GLM-PCA
#'   (default = 1). Ignored if \code{object} is a matrix.
#'
#' @param L the desired number of latent dimensions (integer).
#' @param fam character describing the likelihood to use for the data
#'   (poisson, negative binomial, binomial approximation to multinomial,
#'   bernoulli).
#' @param ctl a list of control parameters for optimization.
#' @param penalty the L2 penalty for the latent factors (default = 1).
#'   Regression coefficients are not penalized.
#' @param verbose logical value indicating whether the current deviance should
#'   be printed after each iteration (default = FALSE).
#' @param init a list containing initial estimates for the factors (\code{U})
#'   and loadings (\code{V}) matrices.
#' @param nb_theta see \code{\link[MASS]{negative.binomial}} (nb_theta -> infty
#'   = Poisson).
#' @param X a matrix of column (observations) covariates. Any column with all
#'   same values (eg. 1 for intercept) will be removed. This is because we force
#'   the intercept and want to avoid collinearity.
#' @param Z a matrix of row (feature) covariates, usually not needed.
#' @param sz numeric vector of size factors to use in place of total counts.
#'
#' @details The basic model is \code{R = AX\'+ZG\'+VU\'}, where \code{E\[Y\] = M
#'   = linkinv(R)}. Regression coefficients are \code{A} and \code{G}, latent
#'   factors are \code{U} and loadings are \code{V}.
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
#' sce <- SingleCellExperiment(assays=list(counts=u))
#' GLMPCA(sce, L = 2)
#' 
#' @import glmpca
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom methods is
#' @export
setMethod(f = "GLMPCA",
		  signature = signature(object = "SummarizedExperiment"),
		  definition = function(object, assay = 1,
		  					  L, fam=c("poi","nb","mult","bern"),
		  					  ctl = list(maxIter=1000, eps=1e-4),
		  					  penalty = 1, verbose = FALSE,
		  					  init = list(factors=NULL, loadings=NULL),
		  					  nb_theta = 1, X = NULL, Z = NULL, sz = NULL){
		  	
		  	res <- glmpca(Y = assay(object, assay),
		  				  L = L, fam = fam, ctl = ctl, penalty = penalty,
		  				  verbose = verbose, init = init, nb_theta = nb_theta,
		  				  X = X, Z = Z, sz = sz)
		  	
		  	if(is(object, 'SingleCellExperiment')){
		  		reducedDim(object, "GLMPCA") <- res$factors
		  	}
		  	object@metadata$glmpca <- res
		  	return(object)
		  })

#' @rdname GLMPCA
#' @importFrom glmpca glmpca
#' @export
setMethod(f = "GLMPCA",
		  signature = signature(object = "matrix"),
		  definition = function(object, assay = 1,
		  					  L, fam=c("poi","nb","mult","bern"),
		  					  ctl = list(maxIter=1000, eps=1e-4),
		  					  penalty = 1, verbose = FALSE,
		  					  init = list(factors=NULL, loadings=NULL),
		  					  nb_theta = 1, X = NULL, Z = NULL, sz = NULL){
		  	
		  	res <- glmpca(Y = object,
		  				  L = L, fam = fam, ctl = ctl, penalty = penalty,
		  				  verbose = verbose, init = init, nb_theta = nb_theta,
		  				  X = X, Z = Z, sz = sz)
		  	
		  	return(res)
		  })
