jeffrey <- function(x, n, conf.level = 0.95){

a <- 1 - conf.level
if (x == 0){ll <- 0}
if (x == n){ul <- 1}
if ((x != 0) & (x != n)){
    ll <- qbeta(a / 2, x + 0.5, n - x + 0.5)  
    ul <- qbeta(1 - a / 2, x + 0.5, n - x + 0.5)}

res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
   
return(res)
}


