confIntIndependentProportion <- function(x, n, conf.level = 0.95)
{
    ## Calculates confidence interval for risk difference of two independent
    ## samples based on individual Wilson intervals

    alpha <- 1 - conf.level

    p1 <- x[1] / n[1]
    p2 <- x[2] / n[2]
    D <- p1 - p2

    ## Newcombe interval
    w1 <- wilson(x[1], n[1], conf.level=conf.level)
    w2 <- wilson(x[2], n[2], conf.level=conf.level)
    D.lower <- D - sqrt((p1 - w1[1]) ^ 2 + (p2 - w2[3]) ^ 2)
    D.upper <- D + sqrt((p1 - w1[3]) ^ 2 + (p2 - w2[1]) ^ 2)
    newcombe <- c(D.lower, D.upper)

    ## Wald interval
    se.D <- sqrt(p1 * (1 - p1) / n[1] + p2 * (1 - p2) / n[2])
    factor <- qnorm(1 - alpha / 2)
    D.lower <- D - factor * se.D
    D.upper <- D + factor * se.D
    wald <- c(D.lower, D.upper)

    res <- list("p1" = p1, "p2" = p2, "d" = as.numeric(D),
                "newcombeCI" = as.numeric(newcombe),
                "waldCI" = as.numeric(wald))
    return(res)
}
