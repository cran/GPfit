
#' Plotting GP model fits
#' 
#' Plots the predicted response and mean squared error (MSE) surfaces for
#' simulators with 1 and 2 dimensional inputs (i.e. d = 1,2).
#' 
#' 
#' @param x a class \code{GP} object estimated by \code{GP_fit}
#' @param M the number of iterations for use in prediction. See
#' \code{\link{predict.GP}}
#' @param range the input range for plotting (default set to \code{[0, 1]})
#' @param resolution the number of points along a coordinate in the specified
#' \code{range}
#' @param colors a vector of length 3 assigning \code{colors[1]} to training
#' design points, \code{colors[2]} to model predictions, and \code{colors[3]}
#' to the error bounds
#' @param line_type a vector of length 2 assigning \code{line_type[1]} to model
#' predictions, and \code{line_type[2]} to the error bounds
#' @param pch a parameter defining the plotting character for the training
#' design points, see `pch' for possible options in \code{\link{par}}
#' @param cex a parameter defining the size of the \code{pch} used for plotting
#' the training design points, see `cex' for possible options in
#' \code{\link{par}}
#' @param legends a parameter that controls the inclusion of a
#' \code{\link{legend}}; by default it is `FALSE'
#' @param surf_check logical, switch between 3d surface and 2d level/contour
#' plotting, the default of \code{FALSE} implies level/contour plotting
#' @param response logical, switch between predicted response and error (MSE)
#' plots, the default of \code{TRUE} displays the response surface
#' @param \dots additional arguments from \code{\link{wireframe}} or
#' \code{\link{levelplot}}
#' @author Blake MacDonald, Hugh Chipman, Pritam Ranjan
#' @seealso \code{\link{GP_fit}} for estimating the parameters of the GP model;
#' \cr \code{\link{predict.GP}} for predicting the response and error surfaces;
#' \cr \code{\link{par}} for additional plotting characters and line types for
#' 1 dimensional plots; \cr \code{\link{wireframe}} and \code{\link{levelplot}}
#' for additional plotting settings in 2 dimensions.
#' @references Ranjan, P., Haynes, R., and Karsten, R. (2011). A
#' Computationally Stable Approach to Gaussian Process Interpolation of
#' Deterministic Computer Simulation Data, Technometrics, 53(4), 366 - 378.
#' @examples
#' 
#' ## 1D Example 1
#' n <- 5
#' d <- 1 
#' computer_simulator <- function(x){
#'     x <- 2 * x + 0.5
#'     y <- sin(10 * pi * x) / (2 * x) + (x - 1)^4
#'     return(y)
#' }
#' set.seed(3)
#' library(lhs)
#' x <- maximinLHS(n,d)
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x,y)
#' plot(GPmodel)
#' 
#' 
#' ## 1D Example 2
#' n <- 7
#' d <- 1
#' computer_simulator <- function(x) {
#'     y <- log(x + 0.1) + sin(5 * pi * x)
#'     return(y)
#' }
#' set.seed(1)
#' library(lhs)
#' x <- maximinLHS(n,d)
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x, y)
#' ## Plotting with changes from the default line type and characters
#' plot(GPmodel, resolution = 100, line_type = c(6,2), pch = 5)
#' 
#' 
#' ## 2D Example: GoldPrice Function
#' computer_simulator <- function(x) {
#'     x1 <- 4 * x[, 1] - 2
#'     x2 <- 4 * x[, 2] - 2
#'     t1 <- 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 
#'         6 * x1 * x2 + 3 * x2^2)
#'     t2 <- 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 
#'         36 * x1 * x2 + 27 * x2^2)
#'     y <- t1 * t2
#'     return(y)
#' }
#' n <- 30 
#' d <- 2
#' set.seed(1)
#' x <- lhs::maximinLHS(n, d) 
#' y <- computer_simulator(x)
#' GPmodel <- GP_fit(x, y)
#' ## Basic level plot
#' plot(GPmodel)
#' ## Adding Contours and increasing the number of levels
#' plot(GPmodel, contour = TRUE, cuts = 50, pretty = TRUE)
#' ## Plotting the Response Surface
#' plot(GPmodel, surf_check = TRUE)
#' ## Plotting the Error Surface with color
#' plot(GPmodel, surf_check = TRUE, response = FALSE, shade = TRUE)
#' 
#' @export
#' @method plot GP
#' @importFrom lattice wireframe levelplot
#' @importFrom graphics legend lines matplot plot

plot.GP <- function(
    x,
    M = 1,
    range = c(0, 1),
    resolution = 50,
    colors = c('black', 'blue', 'red'),
    line_type = c(1, 2),
    pch = 20,
    cex = 1,
    legends = FALSE,
    surf_check = FALSE,
    response = TRUE,
    ...) {
    if (!is.GP(x)) {
        stop("The object in question is not of class \"GP\" \n")
    }
    X <- x$X
    Y <- x$Y
    n <- nrow(X)
    d <- ncol(X)
    
    switch(d,
        {
            # TODO replace with lattice plot
            # to avoid mixed returns
            # make the prediction values
            xvec <- matrix(seq(
                    from = range[1L],
                    to = range[2L],
                    length.out = resolution),
                ncol = 1L)
            GPprediction <- predict.GP(object = x, xnew = xvec, M = M)
            Y_hat <- GPprediction$Y_hat
            MSE <- GPprediction$MSE
            
            # Finding a good range to plot over
            max_height <- max(Y_hat, Y, Y_hat + 2L * sqrt(MSE), Y_hat - 2L * sqrt(MSE))
            min_height <- min(Y_hat, Y, Y_hat + 2L * sqrt(MSE), Y_hat - 2L * sqrt(MSE))
            max_length <- max(X, xvec)
            min_length <- min(X, xvec)
            
            # Plotting the Prediction and errors
            leg.txt <- c(
                expression("Model Prediction:     " ~ hat(y)(x)), 
                expression("Uncertanity Bounds: " ~ hat(y)(x) %+-% 2L %*% s(x)),
                "Design Points")
            matplot(X, Y, 
                cex = cex[1L], col = colors[1L], pch = pch[1L],
                ylim = c(min_height, max_height), 
                xlim = c(min_length, max_length),
                xlab = "x (Input Variable)", 
                ylab = "Model Prediction")
            # Predicted Function
            lines(x = xvec, y = Y_hat, 
                col = colors[2L], lty = line_type[1L])
            # Errors
            lines(x = xvec, y = Y_hat - 2L * sqrt(MSE), 
                col = colors[3L], lty = line_type[2L])
            lines(x = xvec, y = Y_hat + 2L * sqrt(MSE), 
                col = colors[3L], lty = line_type[2L])
            if (legends) {
                legend(x = min_length, y = max_height, 
                    legend = leg.txt, 
                    col = c(colors[2L], colors[3L], colors[1L]), 
                    lty = c(line_type[1L], line_type[2L], -1L), 
                    pch = c(-1L, -1L, pch[1L]), 
                    pt.cex = cex[1L])
            }
        }, 
        {
            # Making the predictions
            xvector <- seq(from = range[1L], to = range[2L], 
                length.out = resolution)
            xvec <- expand.grid(x = xvector, y = xvector)
            xvec <- as.matrix(xvec)
            GPprediction <- predict.GP(object = x, xnew = xvec, M = M)
            Y_hat <- GPprediction$Y_hat
            MSE <- GPprediction$MSE
            dim(Y_hat) <- c(length(xvector), length(xvector))
            dim(MSE) <- c(length(xvector), length(xvector))
            
            if (surf_check) {
                ## Wireframe Plots
                if (response) {
                    h1 <- wireframe(Y_hat,
                        scales = list(arrows = FALSE),
                        row.values = xvector,
                        column.values = xvector,
                        xlab = expression(X[1L]), 
                        ylab = expression(X[2L]),
                        zlab = list('Model Prediction', rot = 90L), ...)
                    print(h1)
                } else {
                    h2 <- wireframe(MSE, 
                        scales = list(arrows = FALSE),
                        row.values = xvector,
                        column.values = xvector,
                        xlab = expression(X[1L]), 
                        ylab = expression(X[2L]),
                        zlab = list('MSE', rot = 90L), ...)
                    print(h2)
                }
            } else {
                ## Contour Plotting
                if (response) {
                    h1 <- levelplot(Y_hat, 
                        xlab = expression(X[1L]), 
                        ylab = expression(X[2L]), 
                        row.values = xvector, 
                        column.values = xvector,
                        xlim = range, ylim = range, ...)
                    print(h1)
                } else {
                    h2 <- levelplot(MSE, 
                        xlab = expression(X[1L]),
                        ylab = expression(X[2L]), 
                        row.values = xvector, 
                        column.values = xvector, 
                        xlim = range, ylim = range, ...)
                    print(h2)
                }
            }
        }, 
        # TODO implement dimension handling
        stop("can not plot in higher than 2 dimensions.\n"))
}
