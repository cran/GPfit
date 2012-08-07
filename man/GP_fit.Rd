\name{GP_fit}
\alias{GP_fit}
\title{Gaussian Process Model fitting}
\description{For an (\code{n x d}) design matrix, \code{X}, and the corresponding
(\code{n x 1}) simulator output \code{Y}, this function fits the GP model and
returns the parameter estimates. The optimization routine assumes that the
inputs are scaled in \eqn{[0,1]^d}{[0,1]^d}.}
\usage{GP_fit(X, Y, control = c(200*d, 80*d, 2*d), nug_thres = 20, 
trace = FALSE, maxit = 100)}
\arguments{
  \item{X}{the (\code{n x d}) design matrix}
  \item{Y}{the (\code{n x 1}) vector of simulator outputs.}
  \item{control}{A vector of parameters used in the search for optimal
 beta. See `Details'.}
  \item{nug_thres}{A parameter used in computing the nugget. See `Details'.}
  \item{trace}{logical, if \code{TRUE}, will provide information on the optim runs}
  \item{maxit}{the maximum number of iterations within \code{\link{optim}}, defaults to 100}
}
\details{
This function fits the following GP model, \eqn{y(x) = \mu + Z(x)}{y(x) = \mu + Z(x)}, \eqn{x \in [0,1]^{d}}{x in [0,1]^d},
where \eqn{Z(x)}{Z(x)} is a GP with \eqn{E(Z(x_i)) = 0}{E(Z(x_i)) = 0}, \eqn{Var(Z(x_i)) = \sigma^2}{Var(Z(x_i)) = \sigma^2}, 
and \eqn{Cov(Z(x_i),Z(x_j)) = \sigma^2R_{ij}}{Cov(Z(x_i),Z(x_j)) = \sigma^2*R_{ij}}. For computational
stability R is replaced with \eqn{R_{\delta_{lb}} = R + \delta_{lb}I}{R_{\delta_{lb}} = R + \delta_{lb}I} 
where \eqn{\delta_{lb}}{\delta_{lb}} is described in Ranjan et al. (2011).\cr

The parameter estimates are obtained by minimizing the deviance using a multi-start gradient based 
search (L-BFGS-B) algorithm. The starting points are selected using the k-means clustering algorithm 
on a large maximin LHD using -2*log(L) as the objective function. The \code{control} parameter determines 
the quality of the starting points of the L-BFGS-B algoritm. \cr

\code{control} is a vector of three tunable parameters used in the deviance optimization algorithm.
The default values correspond to choosing 2*d clusters (using k-means clustering algorithm)
based on 80*d best points (smallest negative 2*log-likelihood, refer to \code{\link{GP_deviance}})
from a 200*d - point random maximin LHD in \eqn{\beta_{0}}{\beta_{0}}.

The \code{nug_thres} parameter is outlined in Ranjan et al. (2011) and is used in finding 
the lower bound on the nugget (\eqn{\delta_{lb}}{\delta_{lb}}).
}

\value{
an object of class \code{GP} that consists of the parameter estimates (\eqn{\hat{\beta}, \delta_{lb}(\hat{\beta}), \hat{\sigma}^2}{\beta_hat, \delta_{lb}(\beta_hat), \sigma_hat^2}), 
the data (\code{X} and \code{Y}), and the correlation structure used
}

\references{Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable 
Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, 
Technometrics, 53(4), 366 - 378.}
\seealso{
\code{\link{plot.GP}} for plotting in 1 and 2 dimensions; \cr 
\code{\link{predict.GP}} for predicting the response and error surfaces; \cr 
\code{\link{optim}} for information on the L-BFGS-B procedure; \cr
\code{\link{GP_deviance}} for computing  the deviance.
}
\keyword{Gaussian Process Modeling}
\author{Blake MacDonald, Hugh Chipman, Pritam Ranjan}
\examples{
## 1D Example 1
n = 5; d = 1; 
computer_simulator <- function(x){
x = 2*x+0.5;
y = sin(10*pi*x)/(2*x) + (x-1)^4;
return(y)
}
set.seed(3);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel)

## 1D Example 2
n = 7; d = 1;
computer_simulator <- function(x) {
y = log(x+0.1)+sin(5*pi*x);
return(y)
}
set.seed(1);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel, digits = 4)

## 2D Example: GoldPrice Function
computer_simulator <- function(x) {
x1=4*x[,1] - 2; x2=4*x[,2] - 2;
t1 = 1 + (x1 + x2 + 1)^2*(19 - 14*x1 + 3*x1^2 - 14*x2 + 
6*x1*x2 + 3*x2^2);
t2 = 30 + (2*x1 -3*x2)^2*(18 - 32*x1 + 12*x1^2 + 48*x2 - 
36*x1*x2 + 27*x2^2);
y = t1*t2;
return(y)
}
n = 30; d = 2;
set.seed(1);
x = maximinLHS(n,d); 
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel)
}  

