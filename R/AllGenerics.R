
#' @title GLM-PCA
#' @rdname GLMPCA
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
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
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
	name = "devianceResiduals",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("devianceResiduals")
	}
)

