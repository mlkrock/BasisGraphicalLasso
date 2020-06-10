#' @title Basis Graphical Lasso
#' @description Estimates the precision matrix from the basis graphical lasso model.
#' @param y Real-valued data matrix of dimension (number of spatial locations) x (number of realizations).
#' @param locs Matrix of real-valued spatial locations of data of dimension (number of spatial locations) x 2.
#' @param lambda Penalty parameter. Can be either a nonnegative real, or a matrix of nonnegative reals whose dimension is the same as the number of graph nodes.
#' @param zero.diagonal.penalty Boolean. If TRUE with a scalar penalty \code{lambda}, then diagonals of the precision matrix are not penalized, Default: TRUE.
#' @param basis Character string for type of basis desired, currently only supports LatticeKrig-type basis.
#' @param Phi Basis matrix, if not specified in \code{basis} option. Rows index location and columns index basis function, Default: NULL.
#' @param distance.penalty If using LatticeKrig Wendland basis functions and a constant lambda, then multiply the distance matrix of the basis centers times lambda for the penalty matrix, Default: FALSE.
#' @param guess An initial guess at the precision matrix. Default: identity matrix.
#' @param tau_sq Nugget variance, estimated by \code{nugget_estimate} if not provided.
#' @param outer_tol Tolerance. Default: see \code{BGL_DC}.
#' @param MAX_ITER Maximum number of iterations. Default: see \code{BGL_DC}.
#' @param MAX_RUNTIME_SECONDS Maximum runtime in seconds. Default: see \code{BGL_DC}.
#' @param verbose Print algorithm details after each iteration. Default: TRUE.
#' @param ... Other options relevant for basis specification such as NC and nlevel for LatticeKrig-type bases.
#' @return Precision matrix of the random coefficients in the weighted sum of basis functions, and estimated (or provided) nugget variance.
#' @details This takes the data itself as input, along with choices or inputs for basis function and nugget variance, and appropriately sends it to the main algorithm in the paper.
#' @examples
#' precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
#'     distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
#'     MAX_RUNTIME_SECONDS=86400, NC=20, nlevel=1)
#'     
#'  # The estimate in paper, with the tolerance set to 1e-2 not 5e-2 and NC=30 not NC=20.
#'  # Takes a few minutes longer than the version above
#'  # precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
#'     # distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
#'     # MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)
#'     
#' # Plot standard errors
#' cholQ <- chol(precision.fit$Q)
#' Phi <- precision.fit$Phi
#' marginal_variances <- rep(NA,dim(Phi)[1])
#' for(i in 1:dim(Phi)[1])
#' {
#'     marginal_variances[i] <- norm(backsolve(cholQ,Phi[i,],transpose=TRUE) ,type="F")^2
#' }
#' quilt.plot(tmin$lon.lat.proj,sqrt(marginal_variances), 
#'      main="Estimated process standard deviation")      
#'      
#' # A (noisy) simulation
#' c_coef <- backsolve(cholQ,rnorm(dim(Phi)[2])) 
#' sim <- Phi %*% c_coef + sqrt(precision.fit$nugget_variance)*rnorm(dim(Phi)[1])
#' quilt.plot(tmin$lon.lat.proj,sim,main="Simulation")     
#' @rdname BGL
#' @export

BGL <- function(y,locs,lambda,zero.diagonal.penalty=TRUE,basis="LatticeKrig",Phi=NULL,guess=NULL,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,tau_sq=NULL,verbose=TRUE,distance.penalty=FALSE,...)
{
  # Set up basis
  basis.setup <- BGLBasisSetup(y=y,locs=locs,basis=basis,Phi=Phi,crossvalidation=FALSE,distance.penalty=distance.penalty,...)
  Phi_Phi <- basis.setup$Phi_Phi
  Phi_S_Phi <- basis.setup$Phi_S_Phi
  trS <- basis.setup$trS
  basiscenters <- basis.setup$basiscenters
  basisdistancematrix <- basis.setup$basisdistancematrix

  # Estimate nugget effect variance
  if(is.null(tau_sq))
  {
    cat('Estimating nugget effect variance: ')
      tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(locs)[1])
      cat(tau_sq, '\n', sep="")
  }

  # Set up identity starting matrix for QUIC algorithm if not prespecified
  l <- dim(Phi_Phi)[1]
  if(is.null(guess))
  {
    guess <- diag(l)
  }

  penaltymatrix <- penaltymatrixsetup(lambda=lambda,zero.diagonal.penalty=TRUE,l=l,basisdistancematrix=basisdistancematrix)

  # Execute DC algorithm
  cat('Running basis graphical lasso algorithm:\n')
  finalguess <- BGL_DC(penaltymatrix,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)

  # Return estimated precision matrix and nugget effect variance
  mylist <- list(Q=finalguess,nugget_variance=tau_sq,Phi=basis.setup$Phi)
  
  return(mylist)
}
