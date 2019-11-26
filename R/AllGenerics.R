
#' @title GLM-PCA
#' @name GLMPCA
#' @export
setGeneric(
	name = "GLMPCA",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("GLMPCA")
	}
)


#' @title Deviance Residuals
#' @name devianceResiduals
#' @export
setGeneric(
	name = "devianceResiduals",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("devianceResiduals")
	}
)

