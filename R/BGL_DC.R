#' @title Basis Graphical Lasso Difference of Convex Algorithm
#' @description Estimates the precision matrix from the basis graphical lasso model.
#' @param lambda Penalty be either a nonnegative real, or a matrix of nonnegative reals whose dimension is the same as the number of graph nodes.
#' @param Phi_Dinv_Phi Inner product of basis matrices, \eqn{\Phi'D^{-1}\Phi}. Typically obtained from using \code{BGLBasisSetup} with \eqn{D = \tau^2 I}.
#' @param Phi_Dinv_S_Dinv_Phi Inner product of the basis matrices and data, \eqn{\Phi'D^{-1}SD^{-1}\Phi}. Typically obtained from using \code{BGLBasisSetup} with \eqn{D = \tau^2 I}. This is where the data directly enters the algorithm. Note: do not compute sample covariance S explicitly. With \code{dat} as the n by m matrix with columns corresponding to realizations of the mean zero spatial field, we have \eqn{S=XX'/m}. So, setting D = I for simplicity, we can compute as \eqn{\Phi'S\Phi} as \code{tcrossprod(crossprod(Phi, dat))/m}.
#' @param guess Initial guess for the precision matrix.
#' @param outer_tol Tolerance. Default: NULL.
#' @param MAX_ITER Maximum number of iterations. Default: NULL.
#' @param MAX_RUNTIME_SECONDS Maximum runtime in seconds. Default: NULL.
#' @param verbose Print algorithm details after each iteration. Default: TRUE.
#' @return Precision matrix of the random coefficients in the weighted sum of basis functions.
#' @details This is the DC algorithm suggested in the paper and is the workhorse of the package. It uses the \code{QUIC} algorithm to solve the inner graphical lasso problem which arises after linearizing all concave functions w.r.t Q in the negative log likelihood. D is the covariance matrix of the additive error term epsilon in the model formulation; in the paper \eqn{D = \tau^2 I}. The package does not estimate non-identity multiple of D but we include the more general version for possible future research.
#' @examples
#' basis.setup <- BGLBasisSetup(y=tmin$data,locs=tmin$lon.lat.proj,basis="LatticeKrig", 
#'      crossvalidation=FALSE,NC=20,nlevel=1)
#' Phi_Phi <- basis.setup$Phi_Phi
#' Phi_S_Phi <- basis.setup$Phi_S_Phi
#' tau_sq <- 2
#' lambda <- matrix(10,nrow=dim(Phi_Phi)[1],ncol=dim(Phi_Phi)[1])
#' diag(lambda) <- 0
#' BGLguess <- BGL_DC(lambda=lambda,Phi_Dinv_Phi=Phi_Phi/tau_sq,
#'      Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2), guess=diag(dim(Phi_Phi)[1]),
#'      outer_tol=0.05,MAX_ITER=50,MAX_RUNTIME_SECONDS=86400)
#' @rdname BGL_DC
#' @export
BGL_DC <- function(lambda,Phi_Dinv_Phi,Phi_Dinv_S_Dinv_Phi,guess,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,verbose=TRUE)
{
  
  if(is.null(outer_tol))
  {
    outer_tol=0.01
  }
  if(is.null(MAX_ITER))
  {
    MAX_ITER=50
  }
  if(is.null(MAX_RUNTIME_SECONDS))
  {
    MAX_RUNTIME_SECONDS=604800
  }
  
  my_norm = 1
  counter = 1L
  start.time <- Sys.time()
  elapsed <- 0
  while (my_norm >= outer_tol && counter <= MAX_ITER && elapsed <= MAX_RUNTIME_SECONDS)
  {
    
    M = chol2inv(chol(guess+Phi_Dinv_Phi))
    S_star <-  M + M %*% Phi_Dinv_S_Dinv_Phi %*% M
    
    new_guess <- QUIC(S_star, rho= lambda,msg=0)$X
    current.time <- Sys.time()
    elapsed <- difftime(current.time,start.time,units="secs")
    my_norm <- norm(new_guess - guess,type="F")/norm(guess,type="F")
    
    if(verbose) {
      cat('Iteration ', counter, '. Relative error: ', my_norm, '. Time elapsed: ', elapsed,"\n", sep="")
    }
    
    guess <- new_guess
    counter <- counter+1
  }
  
  if(my_norm < outer_tol){
    cat("Convergence. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  }  else if(elapsed>MAX_RUNTIME_SECONDS) {
    cat("MAX RUNTIME REACHED. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  } else {
    cat("MAX ITERATIONS REACHED. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  }
  return(guess)
}
