#' Basis Graphical Lasso
#'
#' Estimates a (sparse) precision matrix for the coefficients of the basis functions in the basis graphical lasso model.
#'
#' @docType package
#' @name BasisGraphicalLasso
#' @import QUIC
#' @import LatticeKrig
NULL


#' Minimum temperature residuals
#'
#' The observational minimum temperature residuals dataset from the basis graphical lasso paper. Data matrix is 4577 locations by 150 realizations whose rows are mean-centered; lon.lat are raw longitude/latitude coordinates while lon.lat.proj are sinusoidally-projected coordinates used in the paper's analysis.
#'
#'
#' @docType data
#' @keywords datasets
#' @name tmin
#' @usage data(tmin)
#' @format The observational minimum temperature residuals dataset from the basis graphical lasso paper. Data matrix is 4577 locations by 150 realizations whose rows are mean-centered; lon.lat are raw longitude/latitude coordinates while lon.lat.proj are sinusoidally-projected coordinates used in the paper's analysis.

#' @examples
#' quilt.plot(tmin$lon.lat.proj,tmin$data[,1],main="Day 1",xlab="Easting",ylab="Northing")
"tmin"
