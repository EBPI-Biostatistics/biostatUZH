
confIntRiskDiff <- function(x, n, conf.level = 0.95){

    ci1 <-  wilson(x[1], n[1], conf.level = conf.level)
    ci2 <-  wilson(x[2], n[2], conf.level = conf.level)

    diff <- matrix(ci1[2] - ci2[2])
    se.diff <- sqrt((ci1[2] * (1 - ci1[2])) / n[1] + (ci2[2] * (1 - ci2[2])) / n[2])  
    z <- qnorm((1 + conf.level) / 2)
    wald.lower <- diff - z * se.diff
    wald.upper <- diff + z * se.diff
    score.lower <- diff - sqrt((ci1[2] - ci1[1])^2 + (ci2[3] - ci2[2])^2)
    score.upper <- diff + sqrt((ci2[2] - ci2[1])^2 + (ci1[3] - ci1[2])^2)

    result <- matrix(ncol=2, nrow=2)
    result[,1] <- c(wald.lower, score.lower)
    result[,2] <- c(wald.upper, score.upper)

    out <- data.frame(type = c("Wald", "Wilson"), result)
    names(out) <- c("type", "lower", "upper")

    ret.list <- list("rd" = diff, "CIs" = out)
    return(ret.list)

}


