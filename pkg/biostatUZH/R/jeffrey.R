jeffrey <- function(x, n, conf.level = 0.95)
{
    a <- 1 - conf.level
    
    ll <- if (x == 0) 0 else qbeta(a / 2, x + 0.5, n - x + 0.5)  
    ul <- if (x == n) 1 else qbeta(1 - a / 2, x + 0.5, n - x + 0.5)

    res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
    return(res)
}
