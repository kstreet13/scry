
#' @title GLM-PCA
#' @description This function implements the GLM-PCA dimensionality reduction
#'   method for high-dimensional count data stored in a SingleCellExperiment.
#' @name runGLMPCA
#' 
#' @param object A \link{SingleCellExperiment} object.
#' @param L the desired number of latent dimensions (integer).
#' @param fam character describing the likelihood to use for the data (poisson,
#'   negative binomial, binomial approximation to multinomial, bernoulli).
#' @param ctl a list of control parameters for optimization.
#' @param penalty the L2 penalty for the latent factors (default = 1).
#'   Regression coefficients are not penalized.
#' @param verbose logical value indicating whether the current deviance should
#'   be printed after each iteration (default = FALSE).
#' @param init a list containing initial estimates for the factors (\code{U}) and
#'   loadings (\code{V}) matrices.
#' @param nb_theta see \code{\link[MASS]{negative.binomial}} (nb_theta -> infty
#'   = Poisson).
#' @param X a matrix of column (observations) covariates. Any column with all
#'   same values (eg. 1 for intercept) will be removed. This is because we force
#'   the intercept and want to avoid collinearity.
#' @param Z a matrix of row (feature) covariates, usually not needed.
#' @param sz numeric vector of size factors to use in place of total counts.
#'
#' @details The basic model is \code{R = AX\'+ZG\'+VU\'}, where \code{E\[Y\] = M
#'   = linkinv(R)}. Regression coefficients are \code{A} and \code{G}, latent
#'   factors are \code{U} and loadings are \code{V}.
#' 
#' @return A list containing:
#' 
#' @importFrom glmpca glmpca
#' @export
glmpca <- function(Y, L, fam=c("poi","nb","mult","bern"),
				 ctl = list(maxIter=1000, eps=1e-4),
				 penalty = 1, verbose = FALSE,
				 init = list(factors=NULL, loadings=NULL),
				 nb_theta = 1, X = NULL, Z = NULL, sz = NULL){
  #Y is data with features=rows, observations=cols
  #L is number of desired latent dimensions
  #fam the likelihood for the data 
  #(poisson, negative binomial, binomial approximation to multinomial, bernoulli)
  #ctl a list of control parameters for optimization
  #penalty the L2 penalty for the latent factors
  #regression coefficients are not penalized
  #nb_theta see MASS::negative.binomial (nb_theta->infty = Poisson)
  #X a matrix of column (observations) covariates
  #any column with all same values (eg 1 for intercept) will be removed
  #this is because we force the intercept, so want to avoid collinearity
  #Z a matrix of row (feature) covariates, usually not needed
  #the basic model is R=AX'+ZG'+VU', where E[Y]=M=linkinv(R)
  #regression coefficients are A,G, latent factors are U and loadings V.
  
  fam<-match.arg(fam)
  N<-ncol(Y); J<-nrow(Y)
  #sanity check inputs
  if(fam %in% c("poi","nb","mult","bern")){ stopifnot(min(Y) >= 0) }
  if(fam=="bern"){ stopifnot(max(Y) <= 1) }
  
  #preprocess covariates and set updateable indices
  if(!is.null(X)){ 
    stopifnot(nrow(X)==ncol(Y))
    #we force an intercept, so remove it from X to prevent collinearity
    X<-remove_intercept(X)
    Ko<-ncol(X)+1
  } else {
    Ko<-1
  }
  if(!is.null(Z)){ 
    stopifnot(nrow(Z)==nrow(Y)) 
    Kf<-ncol(Z)
  } else {
    Kf<-0
  }
  lid<-(Ko+Kf)+(1:L)
  uid<-Ko + 1:(Kf+L)
  vid<-c(1:Ko, lid)
  Ku<-length(uid); Kv<-length(vid)
  
  #create glmpca_family object
  gnt<-glmpca_init(Y,fam,sz,nb_theta)
  gf<-gnt$gf; rfunc<-gnt$rfunc; a1<-gnt$intercepts
  
  #initialize U,V, with row-specific intercept terms
  U<-cbind(1, X, matrix(rnorm(N*Ku)*1e-5/Ku,nrow=N))
  if(!is.null(init$factors)){
    #print("initialize factors")
    L0<-min(L,ncol(init$factors))
    U[,(Ko+Kf)+(1:L0)]<-init$factors[,1:L0,drop=FALSE]
  }
  #a1 = naive MLE for gene intercept only
  V<-cbind(a1, matrix(rnorm(J*(Ko-1))*1e-5/Kv,nrow=J))
  V<-cbind(V, Z, matrix(rnorm(J*L)*1e-5/Kv,nrow=J))
  if(!is.null(init$loadings)){
    #print("initialize loadings")
    L0<-min(L,ncol(init$loadings))
    V[,(Ko+Kf)+(1:L0)]<-init$loadings[,1:L0,drop=FALSE]
  }
  
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){
    #rmse[t]<-sd(Y-ilfunc(rfunc(U,V)))
    dev[t]<-gf$dev_func(Y,rfunc(U,V))
    if(t>5 && abs(dev[t]-dev[t-1])/(0.1+abs(dev[t-1]))<ctl$eps){
      break
    }
    if(verbose){ 
      dev_format<-format(dev[t],scientific=TRUE,digits=4)
      print(paste0("Iteration: ",t," | deviance=",dev_format)) 
    }
    #(k %in% lid) ensures no penalty on regression coefficients: 
    for(k in vid){
      ig<- gf$infograd(Y,rfunc(U,V))
      grads<- (ig$grad)%*%U[,k] - penalty*V[,k]*(k %in% lid) 
      infos<- (ig$info) %*% U[,k]^2 + penalty*(k %in% lid)
      V[,k]<-V[,k]+grads/infos
    }
    for(k in uid){
      ig<- gf$infograd(Y,rfunc(U,V))
      grads<- crossprod(ig$grad, V[,k]) - penalty*U[,k]*(k %in% lid) 
      infos<- crossprod(ig$info, V[,k]^2) + penalty*(k %in% lid) 
      U[,k]<-U[,k]+grads/infos
    }
  }
  #postprocessing: include row and columnn labels for regression coefficients
  if(is.null(Z)){
    G<-NULL
  } else {
    G<-U[,Ko+(1:Kf),drop=FALSE]
    rownames(G)<-colnames(Y); colnames(G)<-colnames(Z)
  }
  X<-if(is.null(X)){ matrix(1,nrow=N) } else { cbind(1,X) }
  if(!is.null(colnames(X))){ colnames(X)[1]<-"(Intercept)" }
  A<-V[,1:Ko,drop=FALSE]
  rownames(A)<-rownames(Y); colnames(A)<-colnames(X)
  res<-ortho(U[,lid],V[,lid],A,X=X,G=G,Z=Z,ret="df")
  rownames(res$factors)<-colnames(Y)
  rownames(res$loadings)<-rownames(Y)
  res$dev=dev[1:t]; res$family<-gf
  res
}
