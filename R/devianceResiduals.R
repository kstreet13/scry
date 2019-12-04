
#' @rdname devianceResiduals
#' @description This function computes deviance residuals for count-based
#'   features. Deviance residuals may be used in place of variance for feature
#'   selection.
#' 
#' @param object A \code{\link{SingleCellExperiment}} or
#'   \link{SummarizedExperiment} object. Alternatively,  a matrix of integer
#'   counts.
#' @param assay a character or integer specifying which assay to use for GLM-PCA
#'   (default = 1). Ignored if \code{object} is a matrix.
#' @param fam a character specifying the model type to be used for calculatind
#'   deviance residuals.
#' @param recalcSizeFactors logical, whether column (cell) size factors should
#'   be calculated according to the specified \code{fam} value. If \code{FALSE},
#'   existing size factors are used.
#'   
#' @return The original \code{SingleCellExperiment} or
#'   \code{SummarizedExperiment} object with the deviance residuals, p-values,
#'   and q-values added to the \code{colData}. If the input was a matrix, output
#'   is a \code{data.frame} with these columns.
#'   
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(f = "devianceResiduals",
		  signature = signature(object = "SingleCellExperiment"),
		  definition = function(object, assay = 1,
		  					  fam = c("binomial", "poisson", "geometric"),
		  					  recalcSizeFactors = TRUE){
		  	
		  	fam <- match.arg(fam)
		  	m <- assay(object, assay)
		  	
		  	if(recalcSizeFactors){
		  		sz <- compute_size_factors(m, fam)
		  		sizeFactors(object) <- sz
		  	}else{
		  		sz <- sizeFactors(object)
		  		if(is.null(sz)){
		  			stop('No existing size factors found.')
		  		}
		  	}
		  	
		  	res <- as.data.frame(t(apply(m, 1, .gof_func, sz, fam)))
		  	
		  	res$dev_qval <- p.adjust(res$dev_pval, "BH")
		  	
		  	colData(object) <- cbind(colData(object), res)
		  	
		  	return(object)
		  })	

#' @rdname devianceResiduals
#' @export
setMethod(f = "devianceResiduals",
		  signature = signature(object = "SummarizedExperiment"),
		  definition = function(object, assay = 1,
		  					  fam = c("binomial", "poisson", "geometric"),
		  					  recalcSizeFactors = TRUE){
		  	
		  	fam <- match.arg(fam)
		  	m <- assay(object, assay)
		  	
		  	if(!recalcSizeFactors){
		  		message('Cannot use existing size factors when input is a ',
		  				'SummarizedExperiment. Setting recalcSizeFactors',
		  				'= TRUE.')
		  	}
		  	sz <- compute_size_factors(m, fam)
		  	
		  	res <- as.data.frame(t(apply(m, 1, .gof_func, sz, fam)))
		  	
		  	res$dev_qval <- p.adjust(res$dev_pval, "BH")
		  	
		  	colData(object) <- cbind(colData(object), res)
		  	
		  	return(object)
		  })	

#' @rdname devianceResiduals
#' @export
setMethod(f = "devianceResiduals",
		  signature = signature(object = "matrix"),
		  definition = function(object, assay = 1,
		  					  fam = c("binomial", "poisson", "geometric"),
		  					  recalcSizeFactors = TRUE){
		  	
		  	fam <- match.arg(fam)
		  	m <- object
		  	if(!recalcSizeFactors){
		  		message('Cannot use existing size factors when input is a ',
		  				'matrix. Setting recalcSizeFactors = TRUE.')
		  	}
		  	sz <- compute_size_factors(m, fam)
		  	
		  	res <- as.data.frame(t(apply(m, 1, .gof_func, sz, fam)))
		  	
		  	res$dev_qval <- p.adjust(res$dev_pval, "BH")
		  	
		  	return(res)
		  })	




