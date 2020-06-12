#' @title Compute size factors
#' @name compute_size_factors
#' @description Computes a size factor for each observation (column) of a
#'   count data matrix based on an approximate multinomial model.
#'
#' @param m a matrix or sparse \code{\link{Matrix}} of integer count values.
#' @param fam a string specifying the model type to be used for calculating
#'   size factors. Must be either 'binomial' or 'poisson'.
#'
#' @return A vector of size factors with length equal to the number of columns
#'   of \code{m}.
#'
#' @details Both fam options are approximations to a multinomial model. Size
#'   factors for binomial are simply the column sums (total counts for each
#'   sample). For Poisson, the size factors are given by the column sums after
#'   rescaling to have geometric mean of one. This improves numerical stability.
#'
#' @examples
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' compute_size_factors(u)
#'
#' @importFrom Matrix colSums
#' @importFrom DelayedArray colSums
#' @export
compute_size_factors<-function(m){
    sz<-colSums(m) #base case: binomial
    return(sz)
}
