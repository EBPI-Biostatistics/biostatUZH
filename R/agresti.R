agresti <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)
    
    k <- qnorm(p = (conf.level + 1) / 2)
    
    # Agresti-Coull
    ptilde <- (x + 2) / (n + 4)
    z <- abs(k)
    stderr <- sqrt(ptilde * (1 - ptilde) / (n + 4))
    ll <- max(ptilde - z * stderr, 0)
    ul <- min(ptilde + z * stderr, 1)
    res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
    
    return(res)
}
