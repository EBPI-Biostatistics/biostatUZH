confIntMedian <- function(x, conf.level = 0.95){

    alpha <- 1 - conf.level

    v <- sort(x, na.last = NA)
    n <- length(x)
    m <- median(x)
    
    ## exact using binomial probabilities
    if (n > 0){
        i <- qbinom(alpha / 2, n, 0.5)
        if (i > 0){exact <- c(m, v[i], v[n - i + 1])}
        else {exact <- c(m, NA, NA)}
    } else {exact <- c(NA, NA, NA)
    }
    
    ## approximate according to Altman (2000), p. 37
    approx <- c(NA, NA, NA)
    if (n > 0){
        r <- round(n / 2 - (qnorm(1 - alpha / 2) * sqrt(n) / 2))
        s <- round(1 + n / 2 + (qnorm(1 - alpha / 2) * sqrt(n) / 2))
        approx <- c(m, v[r], v[s])
        }
    
    r <- data.frame(rbind(c(exact[2], exact[1], exact[3]), c(approx[2], approx[1], approx[3])))
    rownames(r) <- c("exact", "approximate")
    colnames(r) <- c("lower", "median", "upper")
    return(r)
}
