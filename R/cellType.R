


#' @rdname predictType
#' @description 
#' 
#' @param target blah 
#' @param reference_dfs 
#' @param types
#' 
#' @return
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
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


