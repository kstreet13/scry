
#' @title GLM-PCA
#' @rdname GLMPCA
#' @export
setGeneric(
	name = "GLMPCA",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("GLMPCA")
	}
)


#' @title Deviance Residuals
#' @rdname devianceResiduals
#' @export
setGeneric(
	name = "devianceResiduals",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("devianceResiduals")
	}
)

