wilson <- function(x, n, conf.level = 0.95){

k <- qnorm(p = (conf.level + 1) / 2)

# Wilson
# function in Hmisc: binconf(9, 20, method = "wilson")                       
res <- (x + k ^ 2 / 2) / (n + k ^ 2) + c(-1, 1) * (k * n ^ 0.5) / (n + k ^ 2) * sqrt(x / n * (1 - x / n) + k ^ 2 / (4 * n))

res <- c("lower" = res[1], "prop" = x / n, "upper" = res[2])
return(res)
}


