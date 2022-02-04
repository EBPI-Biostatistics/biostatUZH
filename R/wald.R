
# Wald confidence interval for a proportion
# similar to Hmisc::binconf(x, n, method="asymptotic")
#
# newly written by LH on 17.02.2016, based on a previous version by KR

wald <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), (x<=n), (n>=1),  conf.level<1, conf.level>0)
    q <- qnorm(p=(1+conf.level)/2)
    pi <- x/n         
    limits <- pi + c(-1, 1)*q*sqrt(pi*(1-pi)/n)
    # NEW: truncate limits at 0 and 1
    res <- c("lower" = max(0, limits[1]), "prop" = pi, "upper" = min(1, limits[2]))
    return(res)
}
