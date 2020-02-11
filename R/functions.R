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

binomial_deviance_residuals<-function(X,p,n){
  #X a matrix, n is vector of length ncol(X)
  stopifnot(length(n)==ncol(X))
  if(length(p)==nrow(X)){ 
    #p is a vector, length must match nrow(X)
    mu<-outer(p,n)
  } else if(!is.null(dim(p)) && dim(p)==dim(X)){
    # p is matrix, must have same dims as X
    mu<-t(t(p)*n)
  } else { stop("dimensions of p and X must match!") }
  term1<-X*log(X/mu)
  term1[is.na(term1)]<-0 #0*log(0)=0
  nx<- t(n-t(X))
  term2<-nx*log(nx/outer(1-p,n))
  #this next line would only matter if all counts
  #were from a single gene, so not checking saves time.
  # term2[is.na(term2)]<-0 
  sign(X-mu)*sqrt(2*(term1+term2))
}

poisson_deviance_residuals<-function(x,xhat){
  #x,xhat assumed to be same dimension
  #sz<-exp(offsets)
  #xhat<-lambda*sz
  term1<-x*log(x/xhat)
  term1[is.na(term1)]<-0 #0*log(0)=0
  s2<-2*(term1-(x-xhat))
  sign(x-xhat)*sqrt(abs(s2))
}

null_residuals<-function(m,fam=c("binomial","poisson"),type=c("deviance","pearson")){
  #m is a matrix
  fam<-match.arg(fam)
  type<-match.arg(type)
  sz<-compute_size_factors(m,fam)
  if(fam=="binomial") {
    phat<-Matrix::rowSums(m)/sum(sz)
    if(type=="deviance"){
      return(binomial_deviance_residuals(m,phat,sz))
    } else { #deviance residuals
      mhat<-outer(phat,sz)
      return((m-mhat)/sqrt(mhat*(1-phat)))
    }
  } else { #fam=="poisson"
    mhat<-outer(Matrix::rowSums(m)/sum(sz), sz) #first argument is "lambda hat" (MLE)
    if(type=="deviance"){ 
      return(poisson_deviance_residuals(m,mhat))
    } else { #pearson residuals
      return((m-mhat)/sqrt(mhat))
    }
  } 
}