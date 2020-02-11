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

poisson_deviance<-function(x,mu,sz){
  #assumes log link and size factor sz on the same scale as x (not logged)
  #stopifnot(all(x>=0 & sz>0))
  2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

binomial_deviance<-function(x,p,n){
  term1<-sum(x*log(x/(n*p)), na.rm=TRUE)
  nx<-n-x
  term2<-sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
  2*(term1+term2)
}

compute_deviance<-function(m,fam=c("binomial","poisson")){
  #m a data matrix with genes=rows
  fam<-match.arg(fam)
  if(is.null(gmeta)){ 
    res<-data.frame()
  } else {
    stopifnot(nrow(m)==nrow(gmeta)) 
  }
  sz<-compute_size_factors(m,fam)
  sz_sum<-sum(sz)
  m<-t(m) #column slicing faster than row slicing for matrix in R.
  #note: genes are now in the COLUMNS of m
  dev_binom<-function(g){
    x<-m[,g]
    binomial_deviance(x,sum(x)/sz_sum,sz)
  }
  dev_poi<-function(g){
    x<-m[,g]
    poisson_deviance(x,sum(x)/sz_sum,sz)
  }
  dev<-if(fam=="binomial"){ dev_binom } else { dev_poi }
  #we can't just use apply() b/c it will coerce a sparse Matrix to dense.
  #note: would be possible to parallelize over cols of m here
  #computation for each gene is independent. Future enhancement?
  vapply(1:ncol(m),dev,FUN.VALUE=0.0) #numeric vector
}
