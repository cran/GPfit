## Iterative Prediction for gaussian process models
##
## May 18th, 2012


#' @name predict
#' @title Model Predictions from GPfit
#' 
#' @description Computes the regularized predicted response \eqn{\hat{y}_{\delta_{lb},M}(x)}
#' and the mean squared error \eqn{s^2_{\delta_{lb},M}(x)} for a new set of
#' inputs using the fitted GP model.
#' 
#' The value of \code{M} determines the number of iterations (or terms) in
#' approximating \eqn{R^{-1} \approx R^{-1}_{\delta_{lb},M}}. The iterative use
#' of the nugget \eqn{\delta_{lb}}, as outlined in Ranjan et al. (2011), is
#' used in calculating \eqn{\hat{y}_{\delta_{lb},M}(x)} and
#' \eqn{s^2_{\delta_{lb},M}(x)}, where \eqn{R_{\delta,M}^{-1} = \sum_{t =
#' 1}^{M} \delta^{t - 1}(R+\delta I)^{-t}}{R_{\delta,M}^{-1} = \sum_{t = 1}^{M}
#' \delta^{t - 1}(R+\delta I)^{-t}}.
#' 
#' @param object a class \code{GP} object estimated by \code{GP_fit}
#' @param xnew the (\code{n_new x d}) design matrix of test points where model
#' predictions and MSEs are desired
#' @param M the number of iterations. See 'Details'
#' @param \dots for compatibility with generic method \code{\link{predict}}
#' @return Returns a list containing the predicted values (\code{Y_hat}), the
#' mean squared errors of the predictions (\code{MSE}), and a matrix
#' (\code{complete_data}) containing \code{xnew}, \code{Y_hat}, and \code{MSE}
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @seealso \code{\link{GP_fit}} for estimating the parameters of the GP model;
#' \cr \code{\link{plot}} for plotting the predicted and error surfaces.
#' @references Ranjan, P., Haynes, R., and Karsten, R. (2011). A
#' Computationally Stable Approach to Gaussian Process Interpolation of
#' Deterministic Computer Simulation Data, Technometrics, 53(4), 366 - 378.
#' @keywords Gaussian Process Model Prediction
#' @importFrom stats predict fitted
#' @examples
#' 
#' ## 1D Example
#' n <- 5
#' d <- 1
#' computer_simulator <- function(x){
#'     x <- 2*x+0.5
#'     sin(10*pi*x)/(2*x) + (x-1)^4
#' }
#' set.seed(3)
#' library(lhs)
#' x <- maximinLHS(n,d)
#' y <- computer_simulator(x)
#' xvec <- seq(from = 0, to = 1, length.out = 10)
#' GPmodel <- GP_fit(x, y)
#' head(fitted(GPmodel))
#' lapply(predict(GPmodel, xvec), head)
#' 
#' 
#' ## 1D Example 2
#' n <- 7
#' d <- 1
#' computer_simulator <- function(x) {
#'     log(x+0.1)+sin(5*pi*x)
#' }
#' set.seed(1)
#' library(lhs)
#' x <- maximinLHS(n,d)
#' y <- computer_simulator(x)
#' xvec <- seq(from = 0,to = 1, length.out = 10)
#' GPmodel <- GP_fit(x, y)
#' head(fitted(GPmodel))
#' predict(GPmodel, xvec)
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
#' n <- 10
#' d <- 2
#' set.seed(1)
#' library(lhs)
#' x <- maximinLHS(n,d) 
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x,y)
#' # fitted values
#' head(fitted(GPmodel))
#' # new data
#' xvector <- seq(from = 0,to = 1, length.out = 10)
#' xdf <- expand.grid(x = xvector, y = xvector)
#' predict(GPmodel, xdf)

NULL

#' @name predict.GP
#' @noRd
#' @export

NULL

#' @describeIn predict The \code{predict} method 
#' returns a list of elements Y_hat (fitted values), 
#' Y (dependent variable), MSE (residuals), and 
#' completed_data (the matrix of independent variables, 
#' Y_hat, and MSE).
#' @export
#' @method predict GP

predict.GP <- function(
    object,
    xnew = object$X,
    M = 1, 
    ...) {
    if (!is.GP(object)) {
        stop("The object in question is not of class \"GP\" \n")
    }
    if (!is.matrix(xnew)) {
        xnew <- as.matrix(xnew)
    }
    if (M <= 0) {
        M <- 1
        warning("M was assigned non-positive, changed to 1. \n")
    }
    
    X <- object$X
    Y <- object$Y
    n <- nrow(X)
    d <- ncol(X)
    corr <- object$correlation_param
    power <- corr$power
    nu <- corr$nu
    
    if (d != ncol(xnew)) {
        stop("The training and prediction data sets are of 
        different dimensions. \n")
    }
    
    beta <- object$beta
    sig2 <- object$sig2
    delta <- object$delta
    
    dim(beta) <- c(d, 1L)
    R <- corr_matrix(X, beta = beta, corr)
    
    Zero <- matrix(0, ncol = 1L, nrow = n)
    One <- matrix(1, ncol = 1L, nrow = n)
    tOne <- matrix(1, ncol = n, nrow = 1L)
    LO <- diag(n)
    Sig <- R + delta * LO
    L <- chol(Sig)
    tL <- t(L)
    
    ## adding in the check on delta to see about the iterative approach
    if (delta == 0) {
        Sig_invOne <- solve(a = L, b = solve(a = tL, b = One))
        Sig_invY <- solve(a = L, b = solve(a = tL, b = Y))
        Sig_invLp <- solve(a = L, b = solve(a = tL, b = LO))
    } else {
    ## Adding in the iterative approach section
        s_Onei <- One
        s_Yi <- Y
        s_Li <- LO
        Sig_invY <- Sig_invOne <- Zero
        Sig_invLp <- matrix(0, ncol = n, nrow = n)
        
        for (it in seq_len(M)) {
            s_Onei <- solve(a = L, b = solve(a = tL, b = delta * s_Onei))
            Sig_invOne <- Sig_invOne + s_Onei/delta
            
            s_Yi <- solve(a = L, b = solve(a = tL, b = delta * s_Yi))
            Sig_invY <- Sig_invY + s_Yi/delta
            
            s_Li <- solve(a = L, b = solve(a = tL, b = delta * s_Li))
            Sig_invLp <- Sig_invLp + s_Li/delta
        }
    }
    
    nnew <- nrow(xnew)
    Y_hat <- rep(0, nnew)
    MSE <- rep(0, nnew)
    ## prediction is different between exponential and 
    ## matern correlation functions
    switch(corr$type,
        "exponential" = {
            for (kk in seq_len(nnew)) {
                ## Changing to accomadate the iterative approach
                xn <- xnew[kk, , drop = FALSE]
                r <- exp(-(abs(X - One %*% xn)^power) %*% (10^beta))
                tr <- t(r)
                yhat <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        tOne + tr) %*% Sig_invY
                Y_hat[kk] <- yhat
                
                ## Adding iterative steps
                if (delta == 0) {
                    Sig_invr <- solve(a = L, b = solve(a = tL, b = r))
                } else {
                    ## if delta != 0, start iterations
                    s_ri <- r
                    Sig_invr <- Zero
                    for (it in seq_len(M)) {
                        s_ri <- solve(a = L, b = solve(a = tL, b = delta * s_ri))
                        Sig_invr <- Sig_invr + s_ri/delta
                    }
                }
                cp_delta_r <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        t(One) + tr) %*% Sig_invr
            
                cp_delta_Lp <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        tOne + tr) %*% Sig_invLp
                mse <- sig2 * (1 - 2 * cp_delta_r + cp_delta_Lp %*% R %*% t(cp_delta_Lp))
                MSE[kk] <- mse * (mse > 0L)
            }
        },
        "matern" = {
            for (kk in seq_len(nnew)) {
                ## Changing to accomodate the iterative approach
                xn <- xnew[kk, , drop = FALSE]
                temp <- 10^beta
                temp <- matrix(temp, 
                    ncol = d, nrow = (length(X) / d), 
                    byrow = TRUE)
                temp <- 2 * sqrt(nu) * abs(X - One %*% xn) * temp
                ID <- which(temp == 0L)
                
                rd <- (1 / (gamma(nu) * 2^(nu - 1))) * (temp^nu) * besselK(temp, nu)    
                rd[ID] <- 1L
                
                r <- matrix(apply(X = rd, MARGIN = 1L, FUN = prod), ncol = 1L)
                tr <- t(r)
                yhat <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        tOne + tr) %*% Sig_invY
                Y_hat[kk] <- yhat
                
                ## Adding iterative steps
                if (delta == 0) {
                    Sig_invr <- solve(a = L, b = solve(a = tL, b = r))
                } else {
                    ## if delta != 0, start iterations
                    s_ri <- r
                    Sig_invr <- Zero
                    for (it in seq_len(M)) {
                        s_ri <- solve(a = L, b = solve(a = tL, b = delta * s_ri))
                        Sig_invr <- Sig_invr + s_ri/delta
                    }
                }
                cp_delta_r <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        tOne + tr) %*% Sig_invr
            
                cp_delta_Lp <- (((1 - tr %*% Sig_invOne) / (tOne %*% Sig_invOne)) %*% 
                        tOne + tr) %*% Sig_invLp
                mse <- sig2 * (1L - 2L * cp_delta_r + cp_delta_Lp %*% R %*% t(cp_delta_Lp))
                MSE[kk] <- mse * (mse > 0L)
            }
        }, 
        stop("unrecognised corr type"))
    
    prediction <- NULL
    
    full_pred <- cbind(xnew, Y_hat, MSE)
    colnames(full_pred) <- c(paste0("xnew.", seq_len(d)), "Y_hat", "MSE")
    prediction$Y_hat <- Y_hat
    prediction$MSE <- MSE
    prediction$complete_data <- full_pred
    
    return(prediction)
}

#' @describeIn predict The \code{fitted} method extracts the complete data.
#' @export
#' @method fitted GP

fitted.GP <- function(object, ...) {
    predict(object)$complete_data
}
