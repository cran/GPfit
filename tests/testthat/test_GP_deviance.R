
context("GP_deviance")

test_that("check GP_deviance", {
    set.seed(6436345)
    nn <- nrow(mtcars)
    x1 <- scale_norm(seq_len(nn))
    y1 <- mtcars$mpg
    gp1 <- GP_deviance(
        beta = 0.9591, 
        X = x1,
        Y = y1)
    expect_equal(
        object = gp1, 
        expected = 255.131524543349)
    
    m2 <- cbind(
        x1,
        mtcars$am, 
        scale_norm(mtcars$disp))
    gp2 <- GP_deviance(
        beta = c(0.9591, 0.1283, 0.8786),
        X = m2,
        Y = y1)
    expect_equal(
        object = gp2, 
        expected = 240.003417343703)
    
    n <- 7
    d <- 1
    computer_simulator <- function(x) {
        log(x + 0.1) + sin(5 * pi * x)
    }
    set.seed(1)
    x <- lhs::maximinLHS(n, d)
    y <- computer_simulator(x)
    beta <- rnorm(1)
    res <- GP_deviance(
        beta = beta,
        X = x,
        Y = y,
        corr = list(
            type = "matern",
            nu = 5/2))
    expect_equal(
        object = res, 
        expected = 37.1900808777125)
    
})
