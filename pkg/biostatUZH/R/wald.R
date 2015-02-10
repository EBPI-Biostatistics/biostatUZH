wald <- function(x, n, conf.level = 0.95){

k <- qnorm(p = (conf.level + 1) / 2)
    
# Wald
# function in Hmisc: binconf(x=9, n=20, method="asymptotic")
pi <- x / n         
res <- pi + c(-1, 1) * k * sqrt(pi * (1 - pi) / n)

res <- c("lower" = res[1], "prop" = x / n, "upper" = res[2])

return(res)
}

