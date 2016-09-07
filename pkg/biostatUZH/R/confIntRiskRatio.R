
confIntRiskRatio <- function(x, n, conf.level = 0.95){
    stopifnot(length(x)==2, length(n)==2, is.wholenumber(x), is.wholenumber(n), (x>0), (x<n),
              conf.level<1, conf.level>0)
    Risk <- x/n
    RiskRatio <- Risk[1]/Risk[2]
    se.log.RiskRatio <- sqrt(sum(1/x)-sum(1/n))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z*se.log.RiskRatio)
    wald.lower <- RiskRatio/EF
    wald.upper <- RiskRatio*EF

    res <- c("lower"=wald.lower, "Risk Ratio"=RiskRatio, "upper"=wald.upper)
    return(res)

}


