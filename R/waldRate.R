
# Wald confidence interval for a Poisson rate
# as in Held, Rufibach, Seifert, Section 8.2
# written by LH on 09.09.2020

#' @export
waldRate <- function(x, t, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), (t>0),  conf.level<1, conf.level>0)
    q <- qnorm(p=(1+conf.level)/2)
    lambda <- x/t
    ef <- exp(q/sqrt(x))
    limits <- c(lambda/ef, lambda*ef)
    res <- c("lower" = limits[1], "rate" = lambda, "upper" = limits[2])
    return(res)
}
