
# calculates limits of Jeffreys confidence interval
# with suitable point estimate
# newly written by LH on 17.02.2016

jeffreys <- function(x, n, conf.level = 0.95){
    q <- (1-conf.level)/2
    alpha <- x + 0.5
    beta <- n - x + 0.5
    point <- qbeta(0.5, alpha, beta)
    ll <- qbeta(q, alpha, beta)
    ul <- qbeta(1 - q, alpha, beta)
    
    res <- c("lower" = ll, "prop" = pi, "upper" = ul)
    return(res)
}
