
# calculates limits of equi-tailed Jeffreys credible interval
# with suitable point estimate (posterior median)
# newly written by LH on 17.02.2016

jeffreys <- function(x, n, conf.level = 0.95){
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x-round(x)) < tol
    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)
    q <- (1-conf.level)/2
    alpha <- x + 0.5
    beta <- n - x + 0.5
    pihat <- qbeta(0.5, alpha, beta)
    limits <- qbeta(c(q, 1-q), alpha, beta)
    res <- c("lower" = limits[1], "pihat" = pihat, "upper" = limits[2])
    return(res)
}
