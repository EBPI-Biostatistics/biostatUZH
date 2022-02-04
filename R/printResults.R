printResults <- function(theta, se.theta, conf.level=0.95, FUN=identity, digits=3){
    z <- qnorm((1+conf.level)/2)
    ci.l <- theta - z*se.theta
    ci.u <- theta + z*se.theta
    p <-  2*pnorm(abs(theta)/se.theta, lower.tail=FALSE)
    effect <- round(FUN(theta), digits=digits)
    ci <- formatCI(c(FUN(ci.l), FUN(ci.u)), digits=digits, text="english")
    p <- formatPval(p)
    res <- c(effect, ci, p)
    names(res) <- list("Effect", 
                       sprintf('%d%% Confidence Interval', conf.level*100),
                       "P-value")
    print(t(res), quote=FALSE)
}
