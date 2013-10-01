##' Fused ANOVA
##'
##' This package is designed to fit accurately the Fused ANOVA model.
##'
##' @section Features:
##'
##' For the moment, only path algorithm with no split is available. 
##'
##' @section Algorithm:
##'
##' Homotopic Algorithm. 
##'
##' @section Technical remarks:
##'
##' Most of the numerical work is done in C++, relying on the
##' \pkg{Rcpp} package. We also provide a cross-validation procedure
##' using the multi-core capability of the computer, through the
##' \pkg{parallel} package. This feature is not available for Windows
##' user, though. 
##'
##' @name fusedanova-package
##' 
##' @docType package
##' @author Pierre Gutierrez \email{pierre.gutierrez@@gmail.com}, Julien chiquet, Guillem Rigaill.
##'
##' @references
##' Shortly coming
##'
##' @import Rcpp parallel ggplot2 grid methods plyr
##' @useDynLib fusedanova
NULL
