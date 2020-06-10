#' @title Evaluate Negative Log Likelihood
#' @description Evaluates the negative log likelihood without the normalizing constant.
#' @param Q Precision matrix.
#' @param Phi_Dinv_Phi Inner product of basis matrices, \eqn{\Phi'D^{-1}\Phi}. Typically obtained from using \code{BGLBasisSetup.R} with \eqn{D = \tau^2 I}.
#' @param Phi_Dinv_S_Dinv_Phi Inner product of the basis matrices and data, \eqn{\Phi'D^{-1}SD^{-1}\Phi}. Typically obtained from using \code{BGLBasisSetup} with \eqn{D = \tau^2 I}. This is where the data directly enters the algorithm. Note: do not compute sample covariance S explicitly. With \code{dat} as the n by m matrix with columns corresponding to realizations of the mean zero spatial field, we have \eqn{S=XX'/m}. So, setting D = I for simplicity, we can compute as \eqn{\Phi'S\Phi} as \code{tcrossprod(crossprod(Phi, dat))/m}.
#' @return -2*log likelihood, ignoring additive constants not depending on Q.
#' @details This returns negative two times the log likelihood, ignoring anything that doesn't depend on Q or \eqn{\tau^2}. Note that the evaluation only requires matrix computations in the number of the basis functions since \eqn{Phi'Phi} and \eqn{Phi'SPhi} are already computed.
#' @rdname unpenalized_likelihood_Q
#' @export

unpenalized_neglikelihood_Q <- function(Q,Phi_Dinv_Phi,Phi_Dinv_S_Dinv_Phi)
{
  cholQ_plus_Phi <- chol(Q + Phi_Dinv_Phi)
  return(2*sum(log(diag(cholQ_plus_Phi))) - 2*sum(log(diag(chol(Q)))) - sum(diag(backsolve(cholQ_plus_Phi,backsolve(cholQ_plus_Phi,Phi_Dinv_S_Dinv_Phi,transpose=TRUE),transpose=FALSE))))
}
