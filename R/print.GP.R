
#' GP model fit Summary
#' 
#' Prints the summary of a class \code{GP} object estimated by \code{GP_fit}
#' 
#' Prints the summary of the class \code{GP} object. It returns the number of
#' observations, input dimension, parameter estimates of the GP model, lower
#' bound on the nugget, and the nugget threshold parameter (described in
#' \code{\link{GP_fit}}).
#' 
#' @param x a class \code{GP} object estimated by \code{GP_fit}
#' @param \dots for compatibility with generic method \code{\link{print}}
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @seealso \code{\link{GP_fit}} for more information on estimating the model;
#' \cr \code{\link{print}} for more description on the \code{print} function.
#' @examples
#' 
#' ## 1D example
#' n <- 5
#' d <- 1 
#' computer_simulator <- function(x){
#'     x <- 2 * x + 0.5
#'     y <- sin(10 * pi * x) / (2 * x) + (x - 1)^4
#'     return(y)
#' }
#' set.seed(3)
#' x <- lhs::maximinLHS(n, d)
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x, y)
#' print(GPmodel)
#' 
#' 
#' ## 2D Example: GoldPrice Function
#' computer_simulator <- function(x) {
#'     x1 <- 4*x[,1] - 2
#'     x2 <- 4*x[,2] - 2
#'     t1 <- 1 + (x1 + x2 + 1)^2*(19 - 14*x1 + 3*x1^2 - 14*x2 + 
#'         6*x1*x2 + 3*x2^2)
#'     t2 <- 30 + (2*x1 -3*x2)^2*(18 - 32*x1 + 12*x1^2 + 48*x2 - 
#'         36*x1*x2 + 27*x2^2)
#'     y <- t1*t2
#'     return(y)
#' }
#' n <- 30 
#' d <- 2
#' set.seed(1)
#' x <- lhs::maximinLHS(n, d) 
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x,y)
#' print(GPmodel, digits = 3)
#' 
#' @export
#' @method print GP

print.GP <- function(x, ...) {
    if (!is.GP(x)) {
        stop("The object in question is not of class \"GP\" \n")
    }
    corr = x$correlation_param
    
    cat("\nNumber Of Observations: n = ")
    cat(nrow(x$X))
    cat("\nInput Dimensions: d = ")
    cat(ncol(x$X))
    cat("\n\n")
    if (corr$type == "exponential"){
        cat("Correlation: Exponential (power = ", corr$power, ")",sep="")
    } else {
        cat("Correlation: Matern (nu = ", corr$nu, ")",sep="")
    }
    cat("\n")
    cat("Correlation Parameters: \n")
    beta_val = data.frame(beta_hat = t(x$beta), row.names = '[1]')
    print(beta_val, ...)
    cat("\n")
    cat("sigma^2_hat: ")
    print(x$sig2, ...)
    cat("\n")
    cat("delta_lb(beta_hat): ")
    new_delta = c(x$delta)
    print(new_delta, ...)
    cat("\n")
    cat("nugget threshold parameter: ")
    cat(x$nugget_threshold_parameter)
    cat("\n\n")
}
