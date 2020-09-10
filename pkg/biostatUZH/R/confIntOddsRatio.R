
confIntOddsRatio <- function(x, n, conf.level = 0.95){
    stopifnot(length(x)==2, length(n)==2, is.wholenumber(x), is.wholenumber(n), (x>0), (x<n),
              conf.level<1, conf.level>0)
    y <- n-x
    Odds <- x/y
    OddsRatio <- Odds[1]/Odds[2]
    se.log.OddsRatio <- sqrt(sum(1/x)+sum(1/y))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z*se.log.OddsRatio)
    wald.lower <- OddsRatio/EF
    wald.upper <- OddsRatio*EF

    res <- c("lower"=wald.lower, "Odds Ratio"=OddsRatio, "upper"=wald.upper)
    return(res)

}


