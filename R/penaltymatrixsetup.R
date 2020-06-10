#' @title Penalty Matrix Setup
#' @description Creates the penalty matrix.
#' @param lambda Scalar or matrix.
#' @param zero.diagonal.penalty Zero penalty along the diagonals of precision. Default: TRUE.
#' @param l Number of basis functions.
#' @param basisdistancematrix The distance between LatticeKrig basis centers. Default: NULL.
#' @return Penalty matrix.
#' @details Returns a penalty matrix.
#' @rdname penaltymatrixsetup
#' @export 
penaltymatrixsetup <- function(lambda,zero.diagonal.penalty=TRUE,l,basisdistancematrix=NULL)
{
  # Set up penalty matrix if not prespecified and do not penalize marginal precision parameters

  if(length(lambda)==1 & !is.null(basisdistancematrix))
  {
    biglambda <- lambda*basisdistancematrix
  }  else if(length(lambda)==1 & is.null(basisdistancematrix))
  {
    biglambda <- matrix(lambda,nrow=l,ncol=l)
  } else {
    biglambda <- lambda
  }
  
  if(zero.diagonal.penalty==TRUE)
  {
    diag(biglambda) <- 0
  }
    
  return(biglambda)
}
