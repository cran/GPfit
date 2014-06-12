\name{corr_matrix}
\alias{corr_matrix}
\title{Power Exponential or Matern Correlation Matrix}
\description{Computes the power exponential or Matern correlation matrix for a set
of \code{n} design points in \code{d}-dimensional input region and a vector of \code{d} 
correlation hyper-parameters.}
\usage{corr_matrix(X, beta, corr=list(type="exponential",power=1.95))}
\arguments{
  \item{X}{the (\code{n x d}) design matrix}
  \item{beta}{a (\code{d x 1}) vector of correlation hyper-parameters in \eqn{(-\infty, \infty)}}
  \item{corr}{a list that specifies the \code{type} of correlation function along with 
the smoothness parameter. The default corresponds to power exponential correlation with 
smoothness parameter "\code{power=1.95}". One can 
instead use \code{corr=list(type = "exponential", power=pk)} with \eqn{1 \le pk \le 2} 
or \code{corr=list(type = "matern", nu=(2*k+1)/2)}, where \eqn{k \in \{0,1,2,...\}}}
}

\value{
The (\code{n x n}) correlation matrix, R, for the design 
matrix (\code{X}) and the hyper-parameters (\code{beta}).
}
\details{
The power exponential correlation function is given by 

\eqn{R_{ij} =  \prod_{k=1}^{d} \exp({-10^{\beta_k}|x_{ik}-x_{jk}|^{power}})}{
R_{ij} = \prod exp(-10^\beta*|x_{i}-x_{j}|^{power})}. \cr


The Matern correlation function given by Santner, Williams, and Notz (2003) is

\eqn{R_{ij} = \prod_{k=1}^{d} \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - 
x_{jk}|\beta_k)^\nu \kappa_{\nu}(2\sqrt{\nu}|x_{ik} - x_{jk}|\beta_k)}{R_{ij} = 
\prod \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - 
x_{jk}|\beta)^\nu \kappa_{\nu}(2\sqrt{\nu}|x_{ik} - x_{jk}|\beta)}, 

where
\eqn{\kappa_{\nu}}{\kappa_{\nu}} is the modified Bessel function of order \eqn{\nu}{\nu}.\cr
}

\references{
MacDonald, K.B., Ranjan, P. and Chipman, H. (2013). GPfit: An R Package for Gaussian 
Process Model Fitting using a New Optimization Algorithm (arXiv:1305.0759).\cr

Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable 
Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, 
Technometrics, 53(4), 366 - 378. \cr

Santner, T.J., Williams, B., and Notz, W. (2003), The design and analysis of 
computer experiments, Springer Verlag, New York.
}

\note{
Both Matern and power exponential correlation functions use the new \eqn{\beta} 
parametrization of hyper-parameters given by \eqn{\theta_k = 10^{\beta_k}} for 
easier likelihood optimization (see MacDonald et al. (2013)).
}

\keyword{Power Exponential Correlation}
\keyword{Gaussian Correlation}
\keyword{Matern Correlation}
\author{Blake MacDonald, Hugh Chipman, Pritam Ranjan}

\examples{
## 1D Example - 1
n = 5; d = 1;
set.seed(3);
library(lhs);
x = maximinLHS(n,d);
beta =  rnorm(1);
R = corr_matrix(x,beta);
print(R)


## 1D Example - 2
beta = rnorm(1)
R = corr_matrix(x,beta,corr = list(type = "matern"))
print(R)


## 2D example - 1
n = 10; d = 2;
set.seed(2);
library(lhs);
x = maximinLHS(n,d); 
beta = rnorm(2);
R = corr_matrix(x,beta,corr = list(type = "exponential", power = 2));
print(R)


## 2D example - 2
beta = rnorm(2)
R = corr_matrix(x,beta,corr = list(type = "matern", nu = 5/2))
print(R)
}
