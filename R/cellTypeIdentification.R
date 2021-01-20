# Get gene-specific distributions
data(params_EM_81020)
pi.all <- params[[1]]
g.on.all <- params[[2]]
g.off.all <- params[[3]]
a.all <- params[[4]]
sigma.all <- params[[5]]
mu <- c(-12.73153,-8.301349)
gaps <- mu[2]+g.on.all-mu[1]-g.off.all
discrim <- rownames(pi.all)[which(gaps>=1&rowSums(pi.all)>=0.05&rowSums(pi.all)<=0.95)]

# Helper function for training reference data
#' @importFrom sads dpoilog
.trainReference <- function(ref) {
  # Subset to genes present in both reference data and knowledge base
  common <- intersect(names(ref),rownames(pi.all))
  N <- sum(ref)
  ref <- ref[common]
  pi.all2 <- pi.all[common,]
  a.all2 <- a.all[common]
  sigma.all2 <- sigma.all[common,]
  g.on.all2 <- g.on.all[common]
  g.off.all2 <- g.off.all[common]
  
  # Compute mixing probability for each component
  prob.exp <- pi.all2[,1]*dnbinom(ref,1,1/((N/(a.all2))+1))
  prob.ln1 <- sapply(common,function(j) 
    pi.all2[j,2]*sads::dpoilog(ref[j],mu[1]+g.off.all2[j]+log(N),sigma.all2[j,1]))
  prob.ln2 <- sapply(common,function(j) 
    (1-rowSums(pi.all2)[j])*sads::dpoilog(ref[j],mu[2]+g.on.all2[j]+log(N),sigma.all2[j,2]))
  return(data.frame(rate=ref/N,exp=prob.exp/(prob.exp+prob.ln1+prob.ln2),ln1=prob.ln1/(prob.exp+prob.ln1+prob.ln2)))
}

#' @title Train reference data to latent states model
#' @rdname trainAllReference
#' @description Fits latent states model to each cell-type present in the provided
#'   reference data by computing cell-type-specific mixing proportions for each
#'   of the three latent states (off-low, off-high, and on). 
#' 
#' @param data an object inheriting from \code{\link{SummarizedExperiment}}
#'   (such as
#'   \code{\link{SingleCellExperiment}}). Alternatively, a matrix or matrix-like
#'   object (such as a sparse \code{\link{Matrix}}) of non-negative integer 
#'   counts, with rownames indicating genes. Genes must be in the format of Ensembl IDs. 
#' @param assay a string or integer specifying which assay contains the count
#'   data (default = 'counts'). Ignored if \code{data} is a matrix-like 
#'   object.
#' @param slot a string or integer specifying which column of \code{colData} 
#'   contains the character cell-type labels for each cell (default = 'labels'). 
#'   Ignored if \code{data} is a matrix-like object.
#' @param labels a vector of characters specifying the cell-type labels for
#'   each cell. Ignored if \code{data} inherits from \code{\link{SummarizedExperiment}}. 
#'
#' @return A list of dataframes indicating the cell-type-specific mixing proportions
#'   for the first two latent states (off-low and off-high) for each gene in
#'   each cell-type.   
#'
#' @details We define three latent states (off-low, indicating unexpressed genes
#'   with very low counts; off-high, indicating unexpressed genes with higher counts;
#'   and on, indicating expression) for each gene in human single-cell data. For
#'   each cell-type in the provided reference data, we learn the cell-type-specific
#'   mixing proportions among these three states for each gene. If there are genes
#'   present in the reference data that are not among the genes in our database,
#'   they are discarded. The object returned here can be passed to 
#'   \code{\link{getBarcode}} to obtain a more interpretable summary of the 
#'   probability that each gene is on in each cell-type, or to \code{\link{classifyTarget}}
#'   to label unknown cells.
#'
#' @examples
#'
#' @references
#' Grabski IN and Irizarry RA (2020). A probabilistic gene barcode for annotation
#' of cell-types from single cell RNA-seq data. \emph{bioRxiv}.
#'
#' @export
setMethod("trainAllReference","matrix",definition = function(data,labels) {
  d.list <- list()
  cell.types <- unique(labels)
  for (c in cell.types) {
    ref <- rowSums(data[,labels==c])
    d.list[[c]] <- .trainReference(ref)
  }
  return(d.list)
})

#' @rdname trainAllReference
#' @importFrom Matrix rowSums
#' @export
setMethod("trainAllReference","Matrix",definition = function(data,labels) {
  d.list <- list()
  cell.types <- unique(labels)
  for (c in cell.types) {
    ref <- rowSums(data[,labels==c])
    d.list[[c]] <- .trainReference(ref)
  }
  return(d.list)
})

#' @rdname trainAllReference
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @export
setMethod("trainAllReference","SummarizedExperiment",definition = function(data,assay="counts",slot="labels") {
  d <- assay(data,assay)
  labels <- colData(data)[,slot]
  cell.types <- unique(labels)
  d.list <- list()
  for (c in cell.types) {
    ref <- rowSums(d[,labels==c])
    d.list[[c]] <- .trainReference(ref)
  }
  return(d.list)
})
 
#' @title Get probabilistic barcode for each cell-type
#' @rdname getBarcode
#' @description Defines probabilistic barcode for each cell-type in reference data
#'   by converting each cell-type's fit to the latent states model into a more
#'   interpretable dataframe of probabilities that each gene is on.
#'
#' @param d.list an object returned from \code{\link{trainAllReference}} 
#'
#' @return A dataframe indicating the probability that each gene is on in each cell-type.
#'
#' @details Fitting each cell-type to the latent states model results in probabilities
#'   of each gene in each cell-type belonging to one of three latent states. This method
#'   takes these fits and computes just the probability that each gene is on (expressed)
#'   in each cell-type, which is a more easily interpretable summary of the model fit. 
#'
#' @examples
#'
#' @references
#' Grabski IN and Irizarry RA (2020). A probabilistic gene barcode for annotation
#' of cell-types from single cell RNA-seq data. \emph{bioRxiv}.
#'
#' @export
setMethod("getBarcode",definition = function(d.list) {
  barcode <- sapply(d.list,function(x) 1-x[,2]-x[,3])
  barcode[barcode<0] <- 0
  colnames(barcode) <- names(d.list)
  rownames(barcode) <- rownames(d.list[[1]])
  return(barcode)
})

#' @title Classify target cells
#' @rdname classifyTarget
#' @description Probabilistically predicts cell-type labels for unknown target cells, 
#'   based on reference data fitted to a latent states model.
#' 
#' @param target an object inheriting from \code{\link{SummarizedExperiment}}
#'   (such as
#'   \code{\link{SingleCellExperiment}}). Alternatively, a matrix or matrix-like
#'   object (such as a sparse \code{\link{Matrix}}) of non-negative integer 
#'   counts, with rownames indicating genes. Genes must be in the format of Ensembl IDs. 
#' @param assay a string or integer specifying which assay contains the count
#'   data (default = 'counts'). Ignored if \code{target} is a matrix-like 
#'   object.
#' @param d.list an object returned from \code{\link{trainAllReference}}
#' @param return.probs if TRUE, the output will be the probabilities that each cell
#'   belongs to each cell-type label under consideration. If FALSE (default), the
#'   output will be the highest probability cell-type label for each cell.
#'
#' @return The inputted \code{\link{SummarizedExperiment}} object with the classification
#'   results appended to the colData. If return.probs = TRUE, the probabilities of each
#'   cell's label will be found under columns named after each cell-type. If 
#'   return.probs = FALSE, the highest probability cell-type label will be found under
#'   a column named predicted_labels.
#'
#' @details We compute the probability of each unknown target cell belonging to each
#'   cell-type label based on the fits of the reference data to a latent states model.
#'   We assume independence of genes, and multiply together the probabilities of each
#'   observed count for each gene in an unknown cell arising under each cell-type label. 
#'   Bayes rule then allows us to compute the probability of each cell-type label. 
#'
#' @examples
#'
#' @references
#' Grabski IN and Irizarry RA (2020). A probabilistic gene barcode for annotation
#' of cell-types from single cell RNA-seq data. \emph{bioRxiv}.
#'
#' @export
setMethod("classifyTarget","matrix",definition = function(target,d.list,return.probs=F) {
   n <- colSums(target)
   return(.computeProbs(target,d.list,n,return.probs))
})

#' @rdname classifyTarget
#' @importFrom Matrix colSums  
#' @export
setMethod("classifyTarget","Matrix",definition = function(target,d.list,return.probs=F) {
   n <- colSums(target)
   return(.computeProbs(target,d.list,n,return.probs))
})

#' @rdname classifyTarget
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<- 
#' @export
setMethod("classifyTarget","SummarizedExperiment",definition = function(target,assay="counts",d.list,return.probs=F) {
  d <- assay(target,assay)
  n <- colSums(d)
  results <- .computeProbs(d,d.list,n,return.probs)
  if (return.probs) {
     colData(target) <- cbind(colData(target),results)
  } else {
     colData(target)[,'predicted_labels'] <- results
  }
  return(target)
}) 

# Compute probabilities for classification of target cells
#' @importFrom sads dpoilog
.computeProbs <- function(target,d.list,n,return.probs=F) {
  genes.s <- intersect(rownames(target),discrim)
  genes.s <- intersect(genes.s,rownames(d.list[[1]]))
  target <- target[genes.s,]
  pi.all2 <- pi.all[genes.s,]
  a2 <- a.all[genes.s]
  sigma.all2 <- sigma.all[genes.s,]
  g.on.all2 <- g.on.all[genes.s]
  g.off.all2 <- g.off.all[genes.s]
  d.list2 <- lapply(seq_len(length(d.list)),function(j) d.list[[j]][genes.s,])
  names(d.list2) <- names(d.list)
  
  # Store probabilities
  gene.components <- array(0,dim=c(length(genes.s),ncol(target),3))
  dimnames(gene.components)[[1]] <- genes.s
  probs <- array(0,dim=c(ncol(target),length(d.list2)))
  
  # For each gene...
  for (j in genes.s) {
    gene.components[j,,1] <- dnbinom(target[j,],1,1/((n/(a2[j]))+1))
    gene.components[j,,2] <- sapply(seq_len(ncol(target)),function(x)
      sads::dpoilog(target[j,x],mu[1]+g.off.all2[j]+log(n[x]),sigma.all2[j,1]))
    gene.components[j,,3] <- sapply(seq_len(ncol(target)),function(x)
      sads::dpoilog(target[j,x],mu[2]+g.on.all2[j]+log(n[x]),sigma.all2[j,2]))
  }
  
  # For each cell-type...
  for (t in seq_len(length(d.list2))) {
    prob.on2 <- 1-rowSums(d.list2[[t]][,2:3])
    prob.on2[prob.on2<0] <- 0
    probs[,t] <- rowSums(log(sweep(t(gene.components[,,1]),2,d.list2[[t]][,2],'*')+
                               sweep(t(gene.components[,,2]),2,d.list2[[t]][,3],'*')+
                               sweep(t(gene.components[,,3]),2,prob.on2,'*')))
  }
  
  if (return.probs) {
    colnames(probs) <- names(d.list2)
    probs <- exp(probs-rowMax(probs))
    return(sweep(probs,1,rowSums(probs),'/'))
  } else {
    return(names(d.list2)[sapply(seq_len(ncol(target)),function(x) which.max(probs[x,]))]) 
  }
}
