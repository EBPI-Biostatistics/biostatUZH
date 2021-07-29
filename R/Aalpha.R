#' @import stats
#' @export
Aalpha <- function(alpha, n, d, MSs, MSe){
    
    ## compute the function A(alpha) on p. 621 of Roussen et al (2003)
    Fa <- qf(alpha, df1 = (n - 1) * (d - 1), df2 = n - 1, ncp = 0)
    res <- (Fa * MSs - MSe) / (d * MSe)
    
    return(res)
}
