#' @title Basis Graphical Lasso Cross Validation
#' @description Estimates the precision matrix from basis graphical lasso model
#' @param y Real-valued data matrix of dimension (number of spatial locations) x (number of realizations).
#' @param kfolds Number of cross validation folds.
#' @param locs Matrix of real-valued spatial locations of data of dimension (number of spatial locations) x 2.
#' @param lambdalist Penalty parameter. Can be either a nonnegative real, or a matrix of nonnegative reals whose dimension is the same as the number of graph nodes.
#' @param zero.diagonal.penalty Boolean, if TRUE with a scalar penalty \code{lambda}, then diagonals of the precision matrix are not penalized. Default: TRUE.
#' @param basis Character string for type of basis desired, currently only supports LatticeKrig-type basis.
#' @param Phi Basis matrix, if not specified in \code{basis} option. Rows index location and columns index basis function, Default: NULL.
#' @param guess An initial guess at the precision matrix. Default of \code{guess} if not specified is the identity matrix.
#' @param distance.penalty If using LatticeKrig Wendland basis functions and a constant lambda, then multiply the distance matrix of the basis centers times lambda for the penalty matrix, Default: FALSE.
#' @param final.guess Return the final estimate with the best penalty parameter, Default: TRUE.
#' @param tau_sq Nugget variance, estimated by \code{nugget_estimate} if not provided.
#' @param outer_tol Tolerance. Default: see \code{BGL_DC}.
#' @param MAX_ITER Maximum number of iterations. Default: see \code{BGL_DC}.
#' @param MAX_RUNTIME_SECONDS Maximum runtime in seconds. Default: see \code{BGL_DC}.
#' @param verbose Print algorithm details after each iteration. Default: TRUE.
#' @param ... Other options relevant for basis specification such as NC and nlevel for LatticeKrig-type bases.
#' @return Precision matrix of the random coefficients in the weighted sum of basis functions.
#' @details This is the algorithm suggested in the paper. It uses the \code{QUIC} algorithm to solve the penalized likelihood problem for each fold, training the model based on "training data" and evaluating fit based on the likelihood with the "testing data". Note that in the paper, the value \eqn{\lambda=7} is reported, which should be consistent with R < version 3.6. With the most recent R version 4.0, the seed is different and \eqn{\lambda=8} is selected.
#' @examples
#' #the full search space considered, commented for runtime, will take very long
#' #CVguess <- BGL_CV(kfolds=2, y=tmin$data, locs=tmin$lon.lat.proj, lambdalist=1:30, basis="LatticeKrig",
#'     #outer_tol=0.05, MAX_ITER=50, MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1,distance.penalty=TRUE)
#' @rdname BGL_CV
#' @export
BGL_CV <- function(kfolds=5,y,locs,lambdalist,zero.diagonal.penalty=TRUE,basis="LatticeKrig",Phi=NULL,guess=NULL,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,tau_sq=NULL,verbose=TRUE,distance.penalty=FALSE,final.guess=TRUE,...)
{
  basis.setup <- BGLBasisSetup(y=y,locs=locs,basis=basis,Phi=Phi,crossvalidation=TRUE,distance.penalty=distance.penalty,...)
  Phi_Phi <- basis.setup$Phi_Phi
  Phi_y <- basis.setup$Phi_y
  Phi_S_Phi <- tcrossprod(Phi_y)/ifelse(ncol(Phi_y)==1,1,ncol(Phi_y)-1)
  trS <- basis.setup$trS
  basiscenters <- basis.setup$basiscenters
  basisdistancematrix <- basis.setup$basisdistancematrix
  
  if(is.null(tau_sq))
  {
    cat('Estimating nugget effect: ')
    tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(locs)[1])
    cat(tau_sq, '\n', sep="")
  }
    
  if(is.null(guess))
  {
    guess <- diag(dim(Phi_Phi)[1])
  }

  set.seed(1)
  folds <- sample(rep(1:kfolds,length.out=ncol(y)))
  likelihood_matrix <- matrix(0,nrow=kfolds,ncol=length(lambdalist))
  
  if(zero.diagonal.penalty==TRUE)
  {  
    zero.diagonal.penalty.list <- vector("list",length=length(lambdalist))
    for(i in 1:length(lambdalist))
    {
      thislambda <- lambdalist[[i]]
      if(length(c(thislambda))==1)
      {
        thislambda <- matrix(thislambda,nrow=dim(Phi_Phi)[1],ncol=dim(Phi_Phi)[1])
        diag(thislambda) <- 0
      } else {
        diag(thislambda) <- 0
      }
      zero.diagonal.penalty.list[[i]] <- thislambda
    }
  }

  cat("Starting CV...")
  
  for(k in 1:kfolds){
     test.indices <- which(folds==k)
     train.indices <- which(folds!=k)

    Phi_Stest_Phi <- tcrossprod(Phi_y[, test.indices])/ ifelse(ncol(Phi_y[,test.indices])==1,1,ncol(Phi_y[,test.indices])-1)
    Phi_Strain_Phi <- tcrossprod(Phi_y[, train.indices])/ ifelse(ncol(Phi_y[,train.indices])==1,1,ncol(Phi_y[,train.indices])-1)

    for(penaltyindex in 1:length(lambdalist))
    {
      penaltymatrix <- penaltymatrixsetup(lambdalist[[penaltyindex]],zero.diagonal.penalty,l,basisdistancematrix)
      foldguess <- BGL_DC(penaltymatrix,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_Strain_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
      likelihood_matrix[k,penaltyindex] <- unpenalized_neglikelihood_Q(foldguess,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_Stest_Phi/(tau_sq^2))
    }
  }
  
  cat("Starting final estimate with best penalty parameter:\n")

  collapse_likelihoods <- colSums(likelihood_matrix)/k
  best_index <- which.min(collapse_likelihoods)
  best_lambda <- penaltymatrixsetup(lambdalist[[best_index]],zero.diagonal.penalty,l,basisdistancematrix)
  
  mylist <- list(nugget_variance=tau_sq,Phi=basis.setup$Phi,best_lambda=best_lambda,best_lambda_index=best_index,likelihood_matrix=likelihood_matrix)
  
  if(final.guess)
  {
    finalguess <- BGL_DC(best_lambda,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
    mylist$Q=finalguess
  }
  
  return(mylist)
  
}
