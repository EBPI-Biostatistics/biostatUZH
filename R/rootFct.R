rootFct <- function(x, alpha, n, d, MSr, MSe){

    Fa <- qf(1 - alpha, df1 = d - 1, df2 = (n - 1) * (d - 1), ncp = x)
    res <- Fa - MSr / MSe
    return(res)
}
