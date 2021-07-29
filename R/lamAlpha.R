#' @export
lamAlpha <- function(alpha, n, d, MSr, MSe){
    
    ## compute the function A(alpha) on p. 621 of Roussen et al (2003)
    limit <- c(10^-5, 1000)
    if (sum(abs(sign(rootFct(limit, alpha, n, d, MSr, MSe)))) == 2){lam.a <- 0} else {
        lam.a <- uniroot(rootFct, interval = limit, alpha, n, d, MSr, MSe)$root}
    return(lam.a)    
}
