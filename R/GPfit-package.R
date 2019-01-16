

#' Gaussian Process Modeling
#' 
#' A computationally stable approach of fitting a Gaussian process (GP) model
#' to simulator outputs. It is assumed that the input variables are continuous
#' and the outputs are obtained from scalar valued deterministic computer
#' simulator.
#' 
#' This package implements a slightly modified version of the regularized GP
#' model proposed in Ranjan et al. (2011). For details see MacDonald et al.
#' (2015). A new parameterization of the Gaussian correlation is used for the
#' ease of optimization. This package uses a multi-start gradient based search
#' algorithm for optimizing the deviance (negative 2*log-likelihood).\cr
#' 
#' For a complete list of functions, use \code{library(help="GPfit")}. \cr The
#' main function for fitting the GP model is \code{\link{GP_fit}}.
#' 
#' @name GPfit-package
#' @aliases GPfit-package GPfit
#' @docType package
#' @author Blake MacDoanld, Hugh Chipman, Pritam Ranjan \cr Maintainer: Hugh
#' Chipman <hugh.chipman@@acadiau.ca>
#' @references MacDonald, K.B., Ranjan, P. and Chipman, H. (2015). GPfit: An R
#' Package for Fitting a Gaussian Process Model to Deterministic Simulator
#' Outputs. Journal of Statistical Software, 64(12), 1-23.
#' \url{http://www.jstatsoft.org/v64/i12/} \cr
#' 
#' Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable
#' Approach to Gaussian Process Interpolation of Deterministic Computer
#' Simulation Data, Technometrics, 53(4), 366 - 378. \cr
#' 
#' Santner, T.J., Williams, B., and Notz, W. (2003), The design and analysis of
#' computer experiments, Springer Verlag, New York. \cr
NULL

#' @title Scale variable into normal range 0, 1
#' @description Perform calculation:
#' (x - min(x)) / (max(x) - min(x))
#' @param x numeric vector
#' @param range numeric vector additional values for shrinking
#' distribution of values within the 0-1 space, without affecting
#' limits of x
#' @return numeric vector
#' @export
#' @examples
#' scale_norm(x = c(-1, 4, 10, 182))
#' # lower bound extended beyond -1
#' # upper bound still range of data
#' scale_norm(x = c(-1, 4, 10, 182), range = c(-100, 100))

scale_norm <- function(x, range = NULL) {
    (x - min(c(x, range), na.rm = TRUE)) / (
        max(c(x, range), na.rm = TRUE) - min(c(x, range), na.rm = TRUE))
}
