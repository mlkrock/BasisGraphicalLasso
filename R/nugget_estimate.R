#' @title Nugget Estimate
#' @description Estimates the nugget variance \eqn{\tau^2} under the simplified assumption that \eqn{Q = \alpha I}. If not desired, user can initially pass \code{tau_sq}.
#' @param Phi_Phi Inner product of basis matrices, \eqn{\Phi'\Phi}. Computed as \code{crossprod(Phi)}, where \code{Phi} has n rows corresponding to locations and l columns corresponding to the basis functions evaluated at those locations.
#' @param Phi_S_Phi Inner product of the basis matrices and data, \eqn{\Phi'S\Phi}. This is where the data directly enters the algorithm. Note: do not compute sample covariance S explicitly. With \code{dat} as the n by m matrix with columns corresponding to realizations of the mean zero spatial field, we have \eqn{S=XX'/m}. So we can compute as \eqn{\Phi'S\Phi} as \code{tcrossprod(crossprod(Phi, dat))/m}.
#' @param trS Trace of empirical covariance matrix.
#' @param n Number of locations.
#' @return Estimate of the nugget variance.
#' @details Performs a joint optimization over \eqn{\alpha,\tau^2} in the -2*negative log likelihood with respect to \eqn{Q = \alpha I} and \eqn{\tau^2}. Calls L-BFBS-B using \code{optim}.
#' @examples
#' basis.setup <- BGLBasisSetup(y=tmin$data,locs=tmin$lon.lat.proj,basis="LatticeKrig",
#'      crossvalidation=FALSE, NC=30,nlevel=1)
#' Phi_Phi <- basis.setup$Phi_Phi
#' Phi_S_Phi <- basis.setup$Phi_S_Phi
#' trS <- basis.setup$trS
#' tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(tmin$lon.lat)[1])
#' @rdname nugget_estimate
#' @export

nugget_estimate <- function(Phi_Phi,Phi_S_Phi,trS,n)
{
  l <- dim(Phi_Phi)[1]

  nugget_function <- function(par)
  {
    diag_parameter <- par[1]
    nug <- par[2]
    Phi_Phi_over_nugget <- Phi_Phi / nug
    Phi_S_Phi_over_nugsq <- Phi_S_Phi / (nug^2)
    cholQ_plus_Phi <- chol(diag_parameter * diag(l) + Phi_Phi_over_nugget)

    return(2*sum(log(diag(cholQ_plus_Phi))) - l*log(diag_parameter) - sum(diag(backsolve(cholQ_plus_Phi,backsolve(cholQ_plus_Phi,Phi_S_Phi_over_nugsq,transpose=TRUE),transpose=FALSE))) + n*log(nug) + trS/nug)
  }

  out <- optim(par=c(1,0.1),fn=nugget_function,hessian=FALSE,method="L-BFGS-B",lower=c(0.01,0.01))
  my_tau_sq <- out$par[2]
  return(my_tau_sq)
}

#same as return(2*sum(log(diag(cholQ_plus_Phi))) - l*log(diag_parameter) - sum(sum(Phi_S_Phi_over_nugsq  * chol2inv(cholQ_plus_Phi))) + n*log(nug) + trS/nug)