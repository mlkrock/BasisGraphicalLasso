#' @title Basis Graphical Lasso Basis Setup
#' @description Sets up basis quadratic forms used in the basis graphical lasso QUIC algorithm.
#' @param y Real-valued data matrix of dimension (number of spatial locations) x (number of realizations); must have at least two realizations.
#' @param locs Matrix of real-valued spatial locations of data of dimension (number of spatial locations) x 2.
#' @param Phi Basis matrix with rows indexing spatial locations and columns indexing basis functions. If not provided, then \eqn{Phi} is generated from the \code{basis} argument. If provided, this calculates the corresponding matrix products.
#' @param basis Character string for type of basis desired, currently only supports LatticeKrig-type basis.
#' @param crossvalidation T/F depending upon whether cross validation is intended; slightly changes what is output.
#' @param distance.penalty T/F will return the locations of the basis centers of the LatticeKrig basis. Then, an example penalty matrix is \code{rdist(basiscenters)}.
#' @param ... Other options relevant for basis specification such as NC and nlevel for LatticeKrig-type bases.
#' @return A list of two matrices \eqn{\Phi'\Phi} and \eqn{\Phi'S\Phi}. Also returns the trace of the sample covariance (for nugget estimation). If cross validation is intended, \eqn{\Phi' Y} is returned instead of \eqn{\Phi' S \Phi}. If 
#' @details Sets up basis inner product matrices that are used in the basis graphical lasso QUIC algorithm.  The first is the inner product of basis matrices, \eqn{\Phi'\Phi}. Computed as \code{crossprod(Phi)}, where \code{Phi} has n rows corresponding to locations and l columns corresponding to the basis functions evaluated at those locations.  Second is inner product of the basis matrices and data, \eqn{\Phi'S\Phi} or \eqn{\Phi'Y} where Y is the data matrix with rows indexing spatial locations. This is where the data directly enters the algorithm. If cross validation is desired, only \eqn{\Phi'Y} is computed and subsequently appropriate columns are used in the cross validation scheme. If \code{distance.penalty=TRUE} then we return the basis centers of the LatticeKrig model for use in the penalty matrix.
#' @examples
#' basis.setup <- BGLBasisSetup(y=tmin$data, locs=tmin$lon.lat.proj, nlevel=1, NC=20)
#' names(basis.setup)
#' @rdname BGLBasisSetup
#' @export

BGLBasisSetup <- function(y,locs,basis="LatticeKrig",Phi=NULL,crossvalidation=FALSE,distance.penalty=FALSE,...)
{
  tmp <- match.call(expand.dots=TRUE)
  if(!is.null(Phi)){
    returnPhi=FALSE
    Phi_Phi <- crossprod(Phi)
    Phi_y <- crossprod(Phi, y)
    trS <- norm(y, "F")^2/ifelse(ncol(y)==1,1,ncol(y)-1)
    
  } else if(is.null(Phi) & basis=="LatticeKrig"){
    returnPhi=TRUE
    if(!all(c("NC","nlevel") %in% names(tmp))){
      stop("Need to specify NC and nlevel for LatticeKrig basis")
    }
    sDomain <- apply(locs, 2, "range")
    # a.wght and nu have no effect on the basis structure but must be specified for LKrigSetup to run
    LKinfo <- LKrigSetup(sDomain, a.wght=4.05, nu=0.5,...)
    Phi <- as.matrix(LKrig.basis(locs, LKinfo))
    Phi_Phi <- crossprod(Phi)
    Phi_y <- crossprod(Phi, y)
    trS <- norm(y, "F")^2/ifelse(ncol(y)==1,1,ncol(y)-1)
  } else {
    stop("LatticeKrig is the only currently supported basis")
  }
  
  if(crossvalidation){
    mylist <- list(Phi_Phi=Phi_Phi, Phi_y=Phi_y, trS=trS)
  } else {
    Phi_S_Phi <- tcrossprod(Phi_y)/ ifelse(ncol(Phi_y)==1,1,ncol(Phi_y)-1)
    mylist <- list(Phi_Phi=Phi_Phi, Phi_S_Phi=Phi_S_Phi, trS=trS)
  }
  
  if(distance.penalty){
    lookvec <- vector('list',LKinfo$nlevel)
    for(i in 1:LKinfo$nlevel){
      look <- LKrigLatticeCenters(LKinfo,Level=i)
      lookvec[[i]] <- expand.grid(look)
    }
    
    basiscenters <- as.matrix(do.call("rbind",lookvec))
    rm(look,lookvec)
    mylist$basiscenters <- basiscenters
    mylist$basisdistancematrix <- rdist(basiscenters)
  }
  
  if(returnPhi) {
    mylist$Phi <- Phi
  }
  
  return(mylist)
}
