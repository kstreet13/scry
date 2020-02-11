#' @title Feature selection by approximate multinomial deviance
#' @rdname devianceFeatureSelection
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
  name = "devianceFeatureSelection",
  signature = 'object',
  def = function(object, ...) {
    standardGeneric("devianceFeatureSelection")
  }
)

#' @title Generalized principal components analysis for non-normally distributed data
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

#' @title Residuals from an approximate multinomial null model
#' @rdname nullResiduals
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
	name = "nullResiduals",
	signature = 'object',
	def = function(object, ...) {
		standardGeneric("nullResiduals")
	}
)
