\name{corr_matrix}
\alias{corr_matrix}
\title{Gaussian Correlation Matrix}
\description{Computes the Gaussian correlation matrix for a set
of design points and a vector of hyper-parameters.}
\usage{corr_matrix(X, beta)}
\arguments{
  \item{X}{the (\code{n x d}) design matrix}
  \item{beta}{a (\code{d x 1}) vector of correlation hyper-parameters}
}
\value{
The (\code{n x n}) correlation matrix, R, for the design matrix (\code{X})
and the hyper-parameters (\code{beta}).
}
\details{
The correlation function is given by 
\eqn{R_{ij} =  \prod_{k=1}^{d} \exp({-10^{\beta_k}(x_{ik}-x_{jk})^2})}{
R_{ij} = \prod exp(-10^\beta*(x_{i}-x_{j})^2)}.
}
\references{Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable 
Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, 
Technometrics, 53(4), 366 - 378.}
\note{Note that a more common formulation would use
\eqn{\theta_{k} = 10^{\beta_{k}}}{\theta_k = 10^\beta_k}, thus \eqn{\beta_{k}}{\beta_k} 
is on a log_10 scale.}
\keyword{Power Exponential Correlation}
\keyword{Gaussian Correlation}
\author{Blake MacDonald, Hugh Chipman, Pritam Ranjan}
\examples{
## 1D Example
n = 5; d = 1;
set.seed(3);
x = maximinLHS(n,d);
beta =  rnorm(1);
R = corr_matrix(x,beta);
print(R)

## 2D example
n = 10; d = 2;
set.seed(1);
x = maximinLHS(n,d); 
beta = rnorm(2);
R = corr_matrix(x,beta);
print(R)
}
