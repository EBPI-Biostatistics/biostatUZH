clopperPearson <- function(x, n, conf.level = 0.95){

    # Compute Clopper-Pearson interval for a binomial proportion
    # Kaspar Rufibach, 2004
    # modified by LH, on 17.02.2016
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x-round(x)) < tol

    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)

    a <- 1 - conf.level
    if (x == 0){
        ll <- 0
        ul <- 1 - (a / 2) ^ (1 / n)}
    else if (x == n){
        ll <- (a / 2) ^ (1 / n)
        ul <- 1
    }
    else {
        ll <- 1/(1 + (n - x + 1) / (x * qf(a / 2, 2 * x, 2 * (n - x + 1))))
        ul <- 1/(1 + (n - x) / ((x + 1) * qf(1 - a / 2, 2 * (x + 1), 2 * (n - x))))
    }
    res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
    return(res)
} 
