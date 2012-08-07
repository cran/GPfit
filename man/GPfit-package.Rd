\name{GPfit-package}
\alias{GPfit-package}
\alias{GPfit}
\docType{package}
\title{
Gaussian Process Modeling
}
\description{
A computationally stable approach of fitting a Gaussian process (GP) 
model to a deterministic simulator.
}
\details{
\tabular{ll}{
Package: \tab GPfit\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2012-08-07\cr
Depends: \tab lhs, lattice\cr
License: \tab GPL-2\cr
}
This package implements a slightly modified version of the GP model proposed 
in Ranjan et al. (2011). A new parameterization of the Gaussian correlation 
is used for the ease of optimization. This package uses a multi-start gradient
based search algorithm for optimizing the deviance (negative 2*log-likelihood).\cr
For a complete list of functions, use \code{library(help="GPfit")}. \cr
The main function for fitting the GP model is \code{\link{GP_fit}}.
}
\author{
Blake MacDoanld, Hugh Chipman, Pritam Ranjan \cr
Maintainer: Hugh Chipman <hugh.chipman@acadiau.ca>
}
\references{
Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable 
Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, 
Technometrics, 53(4), 366 - 378.
}
