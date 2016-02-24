confIntProportion <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)

    res <- data.frame(matrix(NA, ncol = 3))
    colnames(res) <- c("type", "low", "up")

    res[1, 2:3] <- wald(x, n, conf.level = conf.level)[c(1, 3)]
    res[2, 2:3] <- wilson(x, n, conf.level = conf.level)[c(1, 3)]
    res[3, 2:3] <- agresti(x, n, conf.level = conf.level)[c(1, 3)]
    res[4, 2:3] <- jeffreys(x, n, conf.level = conf.level)[c(1, 3)]
    res[5, 2:3] <- clopperPearson(x, n, conf.level = conf.level)[c(1, 3)]

    res[, 1] <- c("Wald", "Wilson", "Agresti", "Jeffreys", "ClopperPearson")

    res <- list("p" = x / n, "CIs" = res)
    return(res)

}
