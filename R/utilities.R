box_cox_inverse <- function(y, lambda = 1)
{
    if (lambda < 0) {
        y[y > -1/lambda] <- NA
    }
    if (lambda == 0) {
        y <- exp(y)
    } else {
        xx <- y * lambda + 1
        y <- sign(xx) * abs(xx)^(1/lambda)
    }
    return(y)
}