
# Wilson confidence interval for a proportion
# similar to Hmisc::binconf(x, n, method = "wilson")
#
# newly written by LH on 17.02.2016, based on a previous version by KR

wilson <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n),
              x<=n, n>=0, # n=0 yields NaN results, but allow for confIntAUC()
              conf.level<1, conf.level>0)
    q <- qnorm(p=(1+conf.level)/2)
    q2 <- q^2
    prop <- x/n
    mid <- (x+q2/2)/(n+q2)
    factor <- (q*sqrt(n))/(n+q2)*sqrt(prop*(1-prop)+q2/(4*n))
    limits <- mid + c(-1,1)*factor
    res <- c("lower"=limits[1], "prop"=prop, "upper"=limits[2])
    return(res)
}
