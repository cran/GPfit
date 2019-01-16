
#' Power Exponential or Matern Correlation Matrix
#' 
#' Computes the power exponential or Matern correlation matrix for a set of
#' \emph{n} design points in \emph{d}-dimensional input region and a vector of
#' \emph{d} correlation hyper-parameters (beta).
#' 
#' The power exponential correlation function is given by
#' 
#' \eqn{R_{ij} = \prod_{k=1}^{d} \exp({-10^{\beta_k}|x_{ik}-x_{jk}|^{power}})}{
#' R_{ij} = \prod exp(-10^\beta*|x_{i}-x_{j}|^{power})}. \cr
#' 
#' The Matern correlation function given by Santner, Williams, and Notz (2003)
#' is
#' 
#' \eqn{R_{ij} = \prod_{k=1}^{d}
#' \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - }{R_{ij} = \prod
#' \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - x_{jk}|10^{\beta_k})^\nu
#' \kappa_{\nu}(2\sqrt{\nu}|x_{ik} -
#' x_{jk}|10^{\beta_k})}\eqn{x_{jk}|10^{\beta_k})^\nu
#' \kappa_{\nu}(2\sqrt{\nu}|x_{ik} - }{R_{ij} = \prod
#' \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - x_{jk}|10^{\beta_k})^\nu
#' \kappa_{\nu}(2\sqrt{\nu}|x_{ik} -
#' x_{jk}|10^{\beta_k})}\eqn{x_{jk}|10^{\beta_k})}{R_{ij} = \prod
#' \frac{1}{\Gamma(\nu)2^{\nu-1}}(2\sqrt{\nu}|x_{ik} - x_{jk}|10^{\beta_k})^\nu
#' \kappa_{\nu}(2\sqrt{\nu}|x_{ik} - x_{jk}|10^{\beta_k})},
#' 
#' where \eqn{\kappa_{\nu}}{\kappa_{\nu}} is the modified Bessel function of
#' order \eqn{\nu}{\nu}.\cr
#' 
#' @param X the (\code{n x d}) design matrix
#' @param beta a (\code{d x 1}) vector of correlation hyper-parameters in
#' \eqn{(-\infty, \infty)}
#' @param corr a list that specifies the \code{type} of correlation function
#' along with the smoothness parameter. The default corresponds to power
#' exponential correlation with smoothness parameter "\code{power=1.95}". One
#' can specify a different power (between 1.0 and 2.0) for the power
#' exponential, or use the Matern correlation function, specified as
#' \code{corr=list(type = "matern", nu=(2*k+1)/2)}, where \eqn{k \in
#' \{0,1,2,...\}}
#' @return The (\code{n x n}) correlation matrix, R, for the design matrix
#' (\code{X}) and the hyper-parameters (\code{beta}).
#' @note Both Matern and power exponential correlation functions use the new
#' \eqn{\beta} parametrization of hyper-parameters given by \eqn{\theta_k =
#' 10^{\beta_k}} for easier likelihood optimization.  That is, \code{beta} is a
#' log scale parameter (see MacDonald et al. (2015)).
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @references MacDonald, K.B., Ranjan, P. and Chipman, H. (2015). 
#' GPfit: An R Package for Fitting a Gaussian Process Model to 
#' Deterministic Simulator Outputs.  
#' Journal of Statistical Software, 64(12), 1-23.
#' \url{http://www.jstatsoft.org/v64/i12/} \cr
#' 
#' Ranjan, P., Haynes, R., and Karsten, R. (2011). 
#' A Computationally Stable
#' Approach to Gaussian Process Interpolation of 
#' Deterministic Computer Simulation Data, 
#' Technometrics, 53(4), 366 - 378. \cr
#' 
#' Santner, T.J., Williams, B., and Notz, W. (2003), 
#' The design and analysis of computer experiments, 
#' Springer Verlag, New York.
#' @keywords Power Exponential Correlation Gaussian Correlation 
#' Matern Correlation
#' @examples
#' ## 1D Example - 1
#' n = 5
#' d = 1
#' set.seed(3)
#' library(lhs)
#' x = maximinLHS(n,d)
#' beta =  rnorm(1)
#' corr_matrix(x,beta)
#' 
#' ## 1D Example - 2
#' beta = rnorm(1)
#' corr_matrix(x,beta,corr = list(type = "matern"))
#' 
#' ## 2D example - 1
#' n = 10
#' d = 2
#' set.seed(2)
#' library(lhs)
#' x = maximinLHS(n,d) 
#' beta = rnorm(2)
#' corr_matrix(x, beta,
#'     corr = list(type = "exponential", power = 2))
#' 
#' ## 2D example - 2
#' beta = rnorm(2)
#' R = corr_matrix(x,beta,corr = list(type = "matern", nu = 5/2))
#' print(R)
#' 
#' @export corr_matrix

corr_matrix <- function(
    X, 
    beta,
    corr = list(
        type = "exponential", 
        power = 1.95)) {
    ## Checking to make sure the data is a matrix, and sets it as one 
    ## if it is not
    if (!is.matrix(X)){
        X = as.matrix(X)
    }
    d = ncol(X)
    n = nrow(X)
    ## Checking the dimensions between the two inputs 
    if (d != length(beta)){
        stop("The dimensions of beta and X do not match. \n")
    }
    # output object
    R <- matrix(0, nrow = n, ncol = n)
    # expand grid
    rcoord <- cbind(
        rep(seq_len(n - 1L), times = rev(seq_len(n - 1L))), 
        unlist(lapply(
            X = rev(seq_len(n - 1L)), 
            FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))

    # X differences
    absdiff <- abs(X[rcoord[, 2L], ] - X[rcoord[, 1L], ])
    # return corr matrix R 
    switch(corr$type, 
        "exponential" = {
            power <- corr$power
            if (is.null(power)){
                power <- 1.95
            }
            
            absdiff <- absdiff^power
            Beta <- matrix(beta, nrow = (length(absdiff) / d), ncol = d, byrow = TRUE)
            Rtemp <- (10^Beta) * absdiff
            # summarise
            Rtemp <- rowSums(Rtemp)
            R[rcoord] <- Rtemp
            R <- R + t(R)
            R <- exp(-R)
            R
        },
        "matern" = {
            nu <- corr$nu
            if (is.null(nu)){
                nu <- 2.5
            }
            
            Beta <- 10^beta
            Beta <- matrix(Beta, ncol = d, nrow = (length(absdiff)/d), byrow = TRUE)
            absdiff <- 2 * sqrt(nu) * absdiff * (Beta)
            pos <- which(absdiff == 0)
            Rtemp <- 1 / (gamma(nu) * 2^(nu - 1)) * (absdiff)^nu * 
                besselK(x = absdiff, nu = nu)
            Rtemp[pos] <- 1
            # summarise
            Rtemp <- apply(X = Rtemp, MARGIN = 1L, FUN = prod)
            # populate matrix upper tri
            R[rcoord] <- Rtemp
            # lower tri
            R <- R + t(R)
            diag(R) <- 1
            R
        }, 
        stop("corr type must be 'exponential' or 'matern'"))
}

