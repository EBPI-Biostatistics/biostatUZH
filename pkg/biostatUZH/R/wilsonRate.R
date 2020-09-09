
# Wilson confidence interval for a rate
# as in Held, Rufibach, Seifert, Section 8.2
# written by LH on 09.09.2020

wilsonRate <- function(x, t, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), t>0, conf.level<1, conf.level>0)
    q <- qnorm(p=(1+conf.level)/2)
    lambda <- x/t
    A <- x + q^2/2
    B <- q/2*sqrt(4*x+q^2)
    limits <- c(A-B, A+B)/t
    res <- c("lower"=limits[1], "rate"=lambda, "upper"=limits[2])
    return(res)
}
