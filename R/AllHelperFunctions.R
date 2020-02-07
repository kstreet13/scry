#' @title compute_size_factors
#' @name compute_size_factors
#' @description This function computes column-wise size factors as appropriate
#'   to a particular, discrete model.
#' 
#' @param m a matrix of integer count values.
#' @param mod a character specifying the model type to be used for calculatind
#'   size factors.
#' 
#' @return A vector of size factors with length equal to the number of columns
#'   of \code{m}.
#' 
#' @importFrom Matrix colSums
#' @export
compute_size_factors <- function(m, 
								 mod = c("binomial", "poisson", "geometric")){
	#given matrix m with samples in the columns
	#compute size factors suitable for the discrete model in 'mod'
	mod <- match.arg(mod)
	sz <- Matrix::colSums(m) #base case, multinomial or binomial
	if(mod == "binomial"){ return(sz) }
	sz <- log(sz)
	sz <- sz - mean(sz) #make geometric mean of sz be 1 for poisson, geometric
	if(mod == "poisson"){ return(exp(sz)) }
	return(sz) #geometric, use log scale size factors
}

#' @importFrom MASS negative.binomial
#' @importFrom stats pchisq glm
.gof_func <- function(x, sz, 
				   mod=c("binomial","poisson","geometric")){
	#Let n=colSums(original matrix where x is a row)
	#if binomial, assumes sz=n, required! So sz>0 for whole vector
	#if poisson, assumes sz=n/geometric_mean(n), so again all of sz>0
	#if geometric, assumes sz=log(n/geometric_mean(n)) which helps numerical stability. Here sz can be <>0
	#note sum(x)/sum(sz) is the (scalar) MLE for "mu" in Poisson and "p" in Binomial
	mod <- match.arg(mod)
	fit <- list(deviance=0,df.residual=length(x)-1,converged=TRUE)
	if(mod=="binomial"){
		fit$deviance <- .binomial_deviance(x,sum(x)/sum(sz),sz)
	} else if(mod=="poisson"){
		fit$deviance <- .poisson_deviance(x,sum(x)/sum(sz),sz)
	} else if(mod=="geometric"){
		if(any(x>0)) {
			fit <- glm(x ~ offset(sz), family=negative.binomial(theta=1))
		}
	} else { stop("invalid model") }
	if(fit$converged){
		dev <- fit$deviance
		df <- fit$df.residual #length(x)-1
		pval <- pchisq(dev, df, lower.tail=FALSE)
		res <- c(dev,pval)
	} else {
		res <- rep(NA,2)
	}
	names(res) <- c("deviance","dev_pval")
	res
}

.poisson_deviance <- function(x, mu, sz){
	#assumes log link and size factor sz on the same scale as x (not logged)
	#stopifnot(all(x>=0 & sz>0))
	2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

.binomial_deviance <- function(x, p, n){
	term1 <- sum(x*log(x/(n*p)), na.rm=TRUE)
	nx <- n-x
	term2 <- sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
	2*(term1 + term2)
}
