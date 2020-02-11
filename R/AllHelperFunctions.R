#' @title Compute size factors
#' @name compute_size_factors
#' @description Computes a size factor for each observation (column) of a 
#'   count data matrix based on an approximate multinomial model.
#' 
#' @param m a matrix or sparse \code{\link{Matrix}} of integer count values.
#' @param fam a string specifying the model type to be used for calculating
#'   size factors.
#' 
#' @return A vector of size factors with length equal to the number of columns
#'   of \code{m}.
#' 
#' @details Both fam options are approximations to a multinomial model. Size factors
#'   for binomial are simply the column sums (total counts for each sample). For
#'   Poisson, the size factors are given by the column sums after rescaling to
#'   have geometric mean of one. This improves numerical stability.
#' 
#' @importFrom Matrix colSums
#' @export
compute_size_factors<-function(m,fam=c("binomial","poisson")){
  #given matrix m with samples in the columns
  #compute size factors suitable for the discrete model in 'fam'
  fam<-match.arg(fam)
  sz<-Matrix::colSums(m) #base case: binomial
  if(fam=="binomial"){ return(sz) }
  #else, Poisson
  lsz<-log(sz)
  #make geometric mean of sz be 1 for poisson
  exp(lsz-mean(lsz))
}

