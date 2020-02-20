


#' @rdname predictType
#' @description This function predicts the identities of each cell in a
#'   single-cell RNAseq dataset, based on reference data corresponding to each
#'   potential cell type of interest.
#'
#' @param target A \code{\link{SummarizedExperiment}} object containing the
#'   single-cell data to be identified.
#' @param reference_dfs A list of reference data fits, returned from the
#'   \code{fitType} function.
#' @param types A character vector indicating the cell type identity of each fit
#'   from \code{reference_dfs}, in order.
#'
#' @return A vector of predicted cell type identities for each cell in
#'   \code{target}. In addition, a new \code{predictedCellType} slot is added to
#'   the \code{SummarizedExperiment} object indicating the cell type identities.
#'   
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom BiocGenerics counts
#' @export
setMethod(f = "predictType",
		  signature = signature(target = "SummarizedExperiment"),
		  definition = function(target, reference_dfs, types){
		  	num_cells <- ncol(target)
		  	num_refs <- length(reference_dfs)
		  	probs <- as.data.frame(
		  		lapply(seq_len(num_refs),function(y){ 
		  			.compute_prob(reference_dfs[[y]], counts(target))}))
		  	chosen <- sapply(seq_len(num_cells),
		  					 function(x){
		  					 	types[which(probs[x,]==max(probs[x,]))]})
		  	colData(target)$predictedCellType <- chosen
		  	return(target)
		  })

#' @rdname fitType
#' @description This function fits a zero-inflated mixture model to reference
#'   single-cell RNA-seq data consisting of entirely one cell type.
#'
#' @param reference A \code{\link{SummarizedExperiment}} object containing the
#'   reference single-cell data.
#' @param iters The number of iterations used in MCMC. Defaults to 5000.
#'
#' @return A \code{data.frame}, where the first column contains the total number
#'   of counts of each gene observed across the entire dataset and the second
#'   column contains the mean parameter for each gene based on the model fit.
#'   
#' @importFrom BiocGenerics counts
#' @export
setMethod(f = "fitType",
		signature = signature(reference = "SummarizedExperiment"),
		definition = function(reference, iters=5000){
			process <- .process_singlecell(counts(reference))
			fit <- .fit_mix(process, iters)
			reference_df <- .prep_df(process, fit)
			return(reference_df)
		})
