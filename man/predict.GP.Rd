\name{predict.GP}
\alias{predict.GP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model Predictions from GPfit}
\description{Computes the predicted response and the mean squared error (MSE) for a set of
new inputs using the fitted GP model.}
\usage{
\method{predict}{GP}(object, xnew = object$X, \dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a class \code{GP} object estimated by \code{GP_fit}}
  \item{xnew}{the (\code{n_new x d}) design matrix of test points where model 
predictions and MSEs are desired}
  \item{\dots}{for compatibility with generic method \code{\link{predict}}}
}
\value{
Returns a list containing the predicted values (\code{Y_hat}), the 
mean squared errors of the predictions (\code{MSE}), and a matrix (\code{complete_data}) containing 
\code{xnew}, \code{Y_hat}, and \code{MSE}}
\references{Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable 
Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, 
Technometrics, 53(4), 366 - 378.}
\keyword{Gaussian Process Model Prediction}
\author{Blake MacDonald, Hugh Chipman, Pritam Ranjan}
\seealso{
\code{\link{GP_fit}} for estimating the parameters of the GP model; \cr
\code{\link{plot.GP}} for plotting the predicted and error surfaces.
}
\examples{
## 1D Example
n = 5; d = 1;
computer_simulator <- function(x){
x = 2*x+0.5;
y = sin(10*pi*x)/(2*x) + (x-1)^4;
return(y)
}
set.seed(3);
x = maximinLHS(n,d);
y = computer_simulator(x);
xvec = seq(from=0,to=1,length.out=10);
GPmodel = GP_fit(x,y);
GPprediction = predict.GP(GPmodel,xvec);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;

## 1D Example 2
n = 7; d = 1;
computer_simulator <- function(x) {
y = log(x+0.1)+sin(5*pi*x);
return(y)
}
set.seed(1)
x = maximinLHS(n,d)
y = computer_simulator(x)
xvec = seq(from=0,to=1,length.out=10);
GPmodel = GP_fit(x,y)
GPprediction = predict.GP(GPmodel,xvec);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;

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
n = 10; d = 2;
set.seed(1);
x = maximinLHS(n,d); 
y = computer_simulator(x);
GPmodel = GP_fit(x,y);

xvector = seq(from=0,to=1,length.out=10);
xvec = expand.grid(x = xvector, y=xvector);
xvec = as.matrix(xvec);
GPprediction = predict.GP(GPmodel,xvec);
yaht = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
}
