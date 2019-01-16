## Optimizes the deviance using multiple optim starts
## 
## May 8th, 2012

#' Gaussian Process Model fitting
#' 
#' For an (\emph{n} x \emph{d}) design matrix, \code{X}, 
#' and the corresponding (\emph{n} x 1) simulator output \code{Y}, 
#' this function fits the GP model and returns the parameter estimates. 
#' The optimization routine assumes that
#' the inputs are scaled to the unit hypercube \eqn{[0,1]^d}.
#' 
#' This function fits the following GP model, 
#' \eqn{y(x) = \mu + Z(x)}{y(x) = \mu + Z(x)}, 
#' \eqn{x \in [0,1]^{d}}{x in [0,1]^d}, where \eqn{Z(x)}{Z(x)} is
#' a GP with mean 0, \eqn{Var(Z(x_i)) = \sigma^2}{Var(Z(x_i)) = \sigma^2}, and
#' \eqn{Cov(Z(x_i),Z(x_j)) = \sigma^2R_{ij}}{Cov(Z(x_i),Z(x_j)) =
#' \sigma^2*R_{ij}}.  Entries in covariance matrix R are determined by
#' \code{corr} and parameterized by \code{beta}, a \code{d}-vector of
#' parameters. For computational stability \eqn{R^{-1}}{R^-1} is replaced with
#' \eqn{R_{\delta_{lb}}^{-1}}, where \eqn{R_{\delta{lb}} = R + \delta_{lb}I}
#' and \eqn{\delta_{lb}} is the nugget parameter described in Ranjan et al.
#' (2011).
#' 
#' The parameter estimate \code{beta} is obtained by minimizing 
#' the deviance using a multi-start gradient based search (L-BFGS-B) 
#' algorithm. The starting points are selected using the k-means 
#' clustering algorithm on a large maximin LHD for values of 
#' \code{beta}, after discarding \code{beta} vectors
#' with high deviance. The \code{control} parameter determines the 
#' quality of the starting points of the L-BFGS-B algoritm.
#' 
#' \code{control} is a vector of three tunable parameters used 
#' in the deviance optimization algorithm. The default values 
#' correspond to choosing 2*d clusters (using k-means clustering 
#' algorithm) based on 80*d best points (smallest deviance, 
#' refer to \code{\link{GP_deviance}}) from a 200*d - point
#' random maximin LHD in \code{beta}. One can change these values 
#' to balance the trade-off between computational cost and robustness 
#' of likelihood optimization (or prediction accuracy).  
#' For details see MacDonald et al. (2015).
#' 
#' The \code{nug_thres} parameter is outlined in Ranjan et al. (2011) and is
#' used in finding the lower bound on the nugget
#' (\eqn{\delta_{lb}}{\delta_{lb}}).
#' 
#' @param X the (\code{n x d}) design matrix
#' @param Y the (\code{n x 1}) vector of simulator outputs.
#' @param control a vector of parameters used in the search for optimal beta
#' (search grid size, percent, number of clusters). See `Details'.
#' @param nug_thres a parameter used in computing the nugget. See `Details'.
#' @param trace logical, if \code{TRUE}, will provide information on the
#' \code{\link{optim}} runs
#' @param maxit the maximum number of iterations within \code{\link{optim}},
#' defaults to 100
#' @param corr a list of parameters for the specifing the correlation to be
#' used. See \code{\link{corr_matrix}}.
#' @param optim_start a matrix of potentially likely starting values for
#' correlation hyperparameters for the \code{\link{optim}} runs, i.e., initial
#' guess of the d-vector \code{beta}
#' @return an object of class \code{GP} containing parameter estimates
#' \code{beta} and \code{sig2}, nugget parameter \code{delta}, the data
#' (\code{X} and \code{Y}), and a specification of the correlation structure
#' used.
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @seealso \code{\link{plot.GP}} for plotting in 1 and 2 dimensions; \cr
#' \code{\link{predict.GP}} for predicting the response and error surfaces; \cr
#' \code{\link{optim}} for information on the L-BFGS-B procedure; \cr
#' \code{\link{GP_deviance}} for computing the deviance.
#' @references MacDonald, K.B., Ranjan, P. and Chipman, H. (2015). GPfit: An R
#' Package for Fitting a Gaussian Process Model to Deterministic Simulator
#' Outputs. Journal of Statistical Software, 64(12), 1-23.
#' \url{http://www.jstatsoft.org/v64/i12/} \cr
#' 
#' Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable
#' Approach to Gaussian Process Interpolation of Deterministic Computer
#' Simulation Data, Technometrics, 53(4), 366 - 378.
#' @keywords Gaussian Process
#' @examples
#' 
#' ## 1D Example 1
#' n = 5
#' d = 1 
#' computer_simulator <- function(x){
#'     x = 2 * x + 0.5
#'     y = sin(10 * pi * x) / (2 * x) + (x - 1)^4
#'     return(y)
#' }
#' set.seed(3)
#' library(lhs)
#' x = maximinLHS(n, d)
#' y = computer_simulator(x)
#' GPmodel = GP_fit(x, y)
#' print(GPmodel)
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
#' GPmodel = GP_fit(x, y)
#' print(GPmodel, digits = 4)
#' 
#' 
#' ## 2D Example: GoldPrice Function
#' computer_simulator <- function(x) {
#'     x1 = 4 * x[, 1] - 2
#'     x2 = 4 * x[, 2] - 2
#'     t1 = 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 
#'         6 * x1 *x2 + 3 * x2^2)
#'     t2 = 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 
#'         36 * x1 * x2 + 27 * x2^2)
#'     y = t1 * t2
#'     return(y)
#' }
#' n = 30
#' d = 2
#' set.seed(1)
#' library(lhs)
#' x = maximinLHS(n, d) 
#' y = computer_simulator(x)
#' GPmodel = GP_fit(x, y)
#' print(GPmodel)
#' 
#' @export GP_fit
#' @importFrom lhs maximinLHS
#' @importFrom stats kmeans optim predict

GP_fit <- function(
    X,
    Y,
    control = c(200 * d, 80 * d, 2 * d),
    nug_thres = 20,
    trace = FALSE,
    maxit = 100,
    corr = list(
        type = "exponential",
        power = 1.95), 
    optim_start = NULL) {
    if (!is.matrix(X)) {
        X = as.matrix(X)
    }
    ## Checking to make sure that the X design is scaled between 0 and 1
    if ((min(X) < 0) | (max(X) > 1) | ((max(X) - min(X)) <= 0.5)) {
        warning("X should be in range (0, 1)")
    }
    n <- nrow(X)
    d <- ncol(X)
    ## Checking the dimensions between the different inputs
    if (n != length(Y)) {
        stop("The dimensions of X and Y do not match. \n")
    }
    if (nug_thres < 10 || nug_thres > 25) {
        warning("nug_thres is outside of the normal range (10, 25)")
    }
    ## Need to make sure that the control values are relatively higher / lower and need 
    ## to print a warning if they are not
    if (length(control) != 3) {
        stop("control is defined incorrectly. Wrong number of arguments. \n")
    }
    if (control[1] < control[2]) { 
        stop("control is defined incorrectly. Need control[1] >= control[2]. \n")
    }
    if (control[2] < control[3]) {
        stop("control is defined incorrectly. Need control[2] >= control[3]. \n")
    }
    ## Checking to see if the vigen starting values  
    ## for the optim runs are in the correct format
    if (!is.null(optim_start)) {
        if (!is.matrix(optim_start)) {
            if (length(optim_start)/d != floor(length(optim_start)/d)) {
                stop("The dimension of optim_start does not match the dimension
                        of the problem \n")
            }
            optim_start = matrix(optim_start, byrow = TRUE, ncol = d)
        } else if (ncol(optim_start) != d) {
            stop("The dimension of optim_start does not match the dimension
                    of the problem \n")
        }
    }
    param_search <- control[1]
    param_percent <- control[2]
    param_clust <- control[3]
    ###########################################################
    ## Using a grid of control[1] points as default, need to find the 
    ## control[2] corresponding values with the lowest deviance
    ## starting out based on the defined range for beta
    beta_range <- switch(corr$type, 
        "exponential" = c((corr$power - 4) - log10(d), log10(5) + corr$power - log10(d)),
        "matern" = c((2) - log10(d), log10(5) + 2 - log10(d)), 
        stop("corr type must be 'exponential' or 'matern'"))
    param_init_ps <- maximinLHS(n = param_search, k = d) * 
        (beta_range[2] - beta_range[1]) + beta_range[1]
    param_lower <- rep(-10, d)
    param_upper <- rep(10, d)
    
    ##-------------------------------------------------------##
    ## Need to evaluate the deviance for the control[1] points
    deviance1 <- cbind(NA, param_init_ps)
    
    deviance1[, 1L] <- apply(
        X = param_init_ps, 
        MARGIN = 1L, 
        FUN = function(row) GP_deviance(
            beta = row,
            X = X,
            Y = Y,
            nug_thres = nug_thres,
            corr = corr))
    
    ## Need to order the initial values based on their deviance
    deviance2 <- deviance1[order(deviance1[, 1L, drop = TRUE]), , drop = FALSE]
    
    ##-------------------------------------------------------##
    ## Taking the control[2] smallest deviance values
    deviance3 <- deviance2[seq_len(param_percent), , drop = FALSE]
    
    ##-------------------------------------------------------##
    ## Going to cluster these control[2] "best" observations into control[3]
    ## groups using K-means, over 5 iterations and the smallest WSS as 
    ## the criterion
    points_percen <- deviance3[, seq_len(d) + 1L, drop = FALSE]
    ## Taking the best of 5 different runs of K-means
    km <- kmeans(
        x = points_percen,
        centers = param_clust, 
        nstart = 5L)
    ##-------------------------------------------------------------------##
    ## Need to set the control[3] best points in cluster as starting values
    param_init <- matrix(0, nrow = param_clust, ncol = d)
    for (ipar in seq_len(param_clust)) {
        id <- which(km$cluster == ipar)
        fid <- which.min(deviance3[id, 1L])
        param_init[ipar, ] <- deviance3[id[fid], seq_len(d) + 1L]
    }
    #############################################################
    ## For the diagonal search: use 3 points along the diagonal, 
    ## 1 near each of the ends and one near the middle of the range, 
    ## this is only necessary above 1 dimension
    if (d >= 2L) {
        param_wrap = matrix(
            c(0.2, 0.5, 0.8) *
                 (beta_range[2L] - beta_range[1L]) + beta_range[1L],
            byrow = TRUE)
    
        # Need to run optim() on the wrapped function the 3 times 
        # in order to find the lowest starting points of the 3 values
        dev <- matrix(NA_real_, nrow = 3L, ncol = 3L)
        for (ipw in seq_len(nrow(param_wrap))) {
            temp <- optim(
                par = param_wrap[ipw],
                fn = dev_wrapper,
                X = X,
                Y = Y,
                nug_thres = nug_thres,
                corr = corr,
                method = "L-BFGS-B",
                lower = param_lower,
                upper = param_upper, 
                control = list(maxit = maxit)) 
            dev[ipw, 1:3] <- c(
                temp$par,
                temp$value, 
                param_wrap[ipw, 1L, drop = TRUE])
        }
        ## Take the best of the 3 based on the likelihood value
        dev <- dev[order(dev[, 2L]), ]
    
    ##---------------------------------------------------------##
        ## Combining the 2*d centers from the clusters with the 
        ## single point from the diagonal search
        param_init <- rbind(rep(dev[1L, 1L], d), param_init)
    }
    
    #############################################################
    ## We now have the control[3]+1 initializing points for the 
    ## optim search, as well as any points that are specified by 
    ## the user
    param_init <- rbind(param_init, optim_start)
    dev_val <- matrix(NA_real_, nrow = nrow(param_init), ncol = d + 1L)
    
    for (ipi in seq_len(nrow(param_init))) {
        temp <- optim(
            par = param_init[ipi, , drop = TRUE],
            fn = GP_deviance,
            X = X,
            Y = Y,
            nug_thres = nug_thres,
            corr = corr,
            method = "L-BFGS-B",
            lower = param_lower,
            upper = param_upper, 
            control = list(maxit = maxit))
        dev_val[ipi, ] <- c(temp$par, temp$value)
    }
    #############################################################
    ## Making a print statement of what the progress of the optimizer
    ## only if trace == TRUE
    if (trace) {
        optim_result = cbind(param_init,dev_val)
        col_name = NULL
        if (d==1){
            row_name = NULL
        } else {
            row_name = c("Diagonal Search")
        }
        for (i in 1:d){
            col_name = cbind(col_name, paste("Beta", as.character(i), "Start"))
        }
        for (i in 1:d){
            col_name = cbind(col_name, paste("Beta", as.character(i), "Final"))
        }
        col_name = cbind(col_name, "Deviance Value")
        for (i in seq_len(param_clust)) {
            row_name = cbind(row_name, paste("Start", as.character(i)))
        }
        colnames(optim_result) = col_name
        rownames(optim_result) = row_name
        print(optim_result)
    }
    # value
    dev_val <- dev_val[order(dev_val[, d + 1L]), ]
    
    beta <- dev_val[1L, seq_len(d), drop = TRUE]
    dim(Y) <- c(n, 1L)
    One <- rep(1, n)
    
    #### sig_invb ####
    
    sigInvb <- sig_invb(X = X, Y = Y, beta = beta, 
        corr = corr, nug_thres = nug_thres)
    
    #### carry on ####
    sig2 <- t(Y - One %*% sigInvb$mu_hat) %*% sigInvb$Sig_invb / n
    
    GP <- list(
        X = X,
        Y = Y,
        sig2 = as.vector(sig2),
        beta = beta,
        delta = sigInvb$delta,
        nugget_threshold_parameter = nug_thres,
        correlation_param = corr)
    
    class(GP) <- "GP"
    return(GP)
}
