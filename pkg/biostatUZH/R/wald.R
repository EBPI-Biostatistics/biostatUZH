
# calculates limits of Wald confidence interval
# newly written by LH on 17.02.2016

wald <- function(x, n, conf.level = 0.95){

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x-round(x)) < tol

    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)
    q <- qnorm(p=(1+conf.level)/2)
    pi <- x/n         
    limits <- pi + c(-1, 1)*q*sqrt(pi*(1-pi)/n)
    res <- c("lower" = max(0, limits[1]), "prop" = pi, "upper" = min(1, limits[2]))
    return(res)
}

