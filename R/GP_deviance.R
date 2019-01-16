
#' Computes the Deviance of a GP model
#' 
#' Evaluates the deviance (negative 2*log-likelihood), as 
#' defined in Ranjan et al. (2011), however the correlation 
#' is reparametrized and can be either power exponential or 
#' Matern as discussed in \code{\link{corr_matrix}}.
#' 
#' 
#' @param beta a (\emph{d} x 1) vector of correlation hyper-parameters, as
#' described in \code{\link{corr_matrix}}
#' @param X the (\emph{n} x \emph{d}) design matrix
#' @param Y the (\emph{n} x 1) vector of simulator outputs
#' @param nug_thres a parameter used in computing the nugget. See
#' \code{\link{GP_fit}}.
#' @param corr a list of parameters for the specifing the correlation to be
#' used. See \code{\link{corr_matrix}}.
#' @return the deviance (negative 2 * log-likelihood)
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @seealso \code{\link{corr_matrix}} for computing the correlation matrix; \cr
#' \code{\link{GP_fit}} for estimating the parameters of the GP model.
#' @references Ranjan, P., Haynes, R., and Karsten, R. (2011). A
#' Computationally Stable Approach to Gaussian Process Interpolation of
#' Deterministic Computer Simulation Data, Technometrics, 53(4), 366 - 378.
#' @keywords Deviance
#' @examples
#' 
#' ## 1D Example 1
#' n = 5
#' d = 1
#' computer_simulator <- function(x) {
#'     x = 2 * x + 0.5
#'     y = sin(10 * pi * x)/(2 * x) + (x - 1)^4
#'     return(y)
#' }
#' set.seed(3)
#' library(lhs)
#' x = maximinLHS(n,d)
#' y = computer_simulator(x)
#' beta =  rnorm(1)
#' GP_deviance(beta,x,y)
#' 
#' 
#' ## 1D Example 2
#' n = 7
#' d = 1
#' computer_simulator <- function(x) {
#'     y <- log(x + 0.1) + sin(5 * pi * x)
#'     return(y)
#' }
#' set.seed(1)
#' library(lhs)
#' x = maximinLHS(n, d)
#' y = computer_simulator(x)
#' beta = rnorm(1)
#' GP_deviance(beta, x, y,
#'     corr = list(type = "matern", nu = 5/2))
#' 
#' ## 2D Example: GoldPrice Function
#' computer_simulator <- function(x) {
#'     x1 = 4 * x[, 1] - 2
#'     x2 = 4 * x[, 2] - 2
#'     t1 = 1 + (x1 + x2 + 1)^2 * 
#'         (19 - 14 * x1 + 3 * x1^2 - 
#'         14 * x2 + 6 * x1 * x2 + 3 * x2^2)
#'     t2 = 30 + (2 * x1 - 3 * x2)^2 * 
#'         (18 - 32 * x1 + 12 * x1^2 + 
#'         48 * x2 - 36 * x1 * x2 + 27 * x2^2)
#'     y = t1 * t2
#'     return(y)
#' }
#' n = 10
#' d = 2
#' set.seed(1)
#' library(lhs)
#' x = maximinLHS(n, d) 
#' y = computer_simulator(x)
#' beta = rnorm(2)
#' GP_deviance(beta, x, y)
#' 
#' @export GP_deviance

GP_deviance <- function(
    beta,
    X,
    Y,
    nug_thres = 20,
    corr = list(
        type = "exponential",
        power = 1.95)) {
    if (!is.matrix(X)) {
        X = as.matrix(X)
    }
    n = nrow(X)
    d = ncol(X)
    ## Checking the dimensions between the different inputs
    if (n != length(Y)) {
        stop("The dimensions of X and Y do not match. \n")
    }
    ## Checking the dimensions of the parameters
    if (d != length(beta)) {
        stop("The dimensions of beta and X do not match \n")
    }
    if (n != length(Y)) {
        stop("The dimensions of X and Y do not match \n")
    }
    
    if (nug_thres < 10 || nug_thres > 25) {
        warning("nug_thres is outside of the normal range (10, 25)")
    }
    
    dim(Y) <- c(n, 1L)
    One <- rep(1, n)
    
    #### sig_invb ####
    
    sigInvb <- sig_invb(X = X, Y = Y, beta = beta, 
        corr = corr, nug_thres = nug_thres)
    
    #### carry on ####
    eig_valL <- eigen(sigInvb$L, only.values = TRUE)$values
    part1 <- 2 * sum(log(abs(eig_valL)))
    # matrix
    part2 <- n * log(t(Y - One %*% sigInvb$mu_hat) %*% sigInvb$Sig_invb)
    devval <- part1 + c(part2)
    
    if (!is.finite(devval)) {
        stop('Infinite values of the Deviance Function, 
            unable to find optimum parameters \n')
    }
    return(devval)
}

#' @title Internal tools
#' @description shared utilities between \code{\link{GP_deviance}}
#' and \code{\link{GP_fit}}
#' @inheritParams GP_deviance
#' @return list with elements delta, L, mu_hat, Sig_invb
#' @examples 
#' set.seed(3234)
#' GPfit:::sig_invb(
#'     X = matrix((0:10) / 10), 
#'     Y = runif(11), 
#'     beta = 1.23)

sig_invb <- function(
    X, Y, beta, 
    corr = list(
        type = "exponential",
        power = 1.95), 
    nug_thres = 20) {
    n <- nrow(X)
    One <- rep(1L, n)
    R <- corr_matrix(X = X, beta = beta, corr = corr)
    temp <- eigen(x = R, symmetric = TRUE, only.values = TRUE)
    eig_val <- temp$values
    condnum <- kappa(z = R, triangular = TRUE, exact = TRUE)
    delta <- max(c(0, abs(eig_val[1L]) * (condnum - exp(nug_thres)) / 
            (condnum * (exp(nug_thres) - 1))), na.rm = TRUE)
    
    LO <- diag(n)
    Sig <- R + delta * LO
    
    L <- chol(x = Sig)
    
    Sig_invOne <- solve(L, solve(t(L), One))
    Sig_invY <- solve(L, solve(t(L), Y))
    
    mu_hat <- solve(t(One) %*% Sig_invOne, t(One) %*% Sig_invY)
    Sig_invb <- solve(L, solve(t(L), (Y - One %*% mu_hat)))
    return(list(
        delta = delta,
        L = L,
        mu_hat = mu_hat, 
        Sig_invb = Sig_invb))
}
