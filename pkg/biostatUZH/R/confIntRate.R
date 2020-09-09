confIntRate <- function(x, t, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), (t>0),  conf.level<1, conf.level>0)

    res <- data.frame(matrix(NA, ncol = 3))
    colnames(res) <- c("type", "lower", "upper")

    res[1, 2:3] <- waldRate(x, t, conf.level = conf.level)[c(1, 3)]
    res[2, 2:3] <- wilsonRate(x, t, conf.level = conf.level)[c(1, 3)]

    res[, 1] <- c("Wald", "Wilson")

    res <- list("rate" = x / t, "CIs" = res)
    return(res)

}
