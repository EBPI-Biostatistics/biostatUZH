confIntPairedProportion <- function(x, conf.level = 0.95){

    ## confidence interval for a paired proportion, according to Newcombe (1998)
    ## see Altman (2000), "Statistics with confidence", p. 52 for a worked out example
    r <- x[1, 1]
    s <- x[1, 2]
    t <- x[2, 1]
    u <- x[2, 2]
    n <- r + s + t + u
    
    p1 <- (r + s) / n
    p2 <- (r + t) / n
    D <- p1 - p2
    
    ## compute Wilson CIs for p1 and p2
    ci1 <- wilson(x = r + s, n = n, conf.level = conf.level)
    ci2 <- wilson(x = r + t, n = n, conf.level = conf.level)
    l1 <- ci1[1]
    u1 <- ci1[3]
    l2 <- ci2[1]
    u2 <- ci2[3]

    phi <- 0
    if (max(r + s, t + u, r + t, s + u) > 0){
        A <- (r + s) * (t + u) * (r + t) * (s + u)
        B <- r * u - s * t
        C <- 0
        if (B > n / 2){C <- B - n / 2}
        if (B < 0){C <- B}
        phi <- C / sqrt(A)    
    }
    
    newcombe1 <- D - sqrt((p1 - l1) ^ 2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2) ^ 2)
    newcombe2 <- D + sqrt((p2 - l2) ^ 2 - 2 * phi * (p2 - l2) * (u1 - p1) + (u1 - p1) ^ 2)
    
    newcombe <- c(newcombe1, newcombe2)
    return(list("p1" = p1, "p2" = p2, "newcombeCI" = as.numeric(newcombe)))
}

