#' @title Feature selection by approximate multinomial deviance
#' @rdname devianceFeatureSelection
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
    name = "devianceFeatureSelection", 
    def = function(object, ...) {
        standardGeneric("devianceFeatureSelection")
    }
)

#' @title Generalized principal components analysis for non-normally distributed
#'   data
#' @rdname GLMPCA
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
    name = "GLMPCA",
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

#' @title Train reference data to latent states model
#' @rdname trainAllReference
#' @param ... for the generic, additional arguments pass to object-specific methods.
#' @export 
setGeneric(name="trainAllReference",def=function(data,...) {
  standardGeneric("trainAllReference")
})

#' @title Get probabilistic barcode for each cell-type
#' @rdname getBarcode
#' @export 
setGeneric(name="getBarcode",def=function(d.list) {
  standardGeneric("getBarcode")
})

#' @title Classify target cells
#' @rdname classifyTarget
#' @param ... for the generic, additional arguments pass to object-specific methods.
#' @export 
setGeneric(name="classifyTarget",def=function(target,...) {
  standardGeneric("classifyTarget")
})

