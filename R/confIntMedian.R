#' Compute exact confidence interval for median of a sample based on order
#' statistics
#' 
#' Given a sample \eqn{X_1, ..., X_n}, this function computes an exact
#' confidence interval for the median of the sample based on the binomial
#' distribution.
#' 
#' 
#' @param x Vector of observations.
#' @param conf.level Confidence level for confidence interval.
#' @return Table containing the median and the confidence interval.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Lehmann, E. (1975). Nonparametrics: Statistical Methods Based on
#' Ranks. \emph{Holden-Day}, 182--183.
#' 
#' A worked out example (in Section 5, p. 38) can be found in
#' 
#' Altman, D.G., Machin, D., Bryant, T.N., Gardner, M.J.(2000). Statistics with
#' confidence. University Press Belfast.
#' @keywords nonparametric
#' @examples
#' 
#' # generate random sample
#' set.seed(1977)
#' n <- 40
#' x <- rnorm(n)
#' 
#' # compute CI for median
#' confIntMedian(x, conf.level = 0.95)
#' 
#' # data from Altman (2000), p. 38
#' x <- c(66, 71.2, 83, 83.6, 101, 107.6, 122, 143, 160, 177, 414)
#' confIntMedian(x, conf.level = 0.95)
#' 
#' @export
confIntMedian <- function(x, conf.level = 0.95){

    alpha <- 1 - conf.level

    v <- sort(x, na.last = NA)
    n <- length(x)
    m <- median(x)
    
    ## exact using binomial probabilities
    if (n > 0){
        i <- qbinom(alpha / 2, n, 0.5)
        if (i > 0){exact <- c(m, v[i], v[n - i + 1])}
        else {exact <- c(m, NA, NA)}
    } else {exact <- c(NA, NA, NA)
    }
    
    ## approximate according to Altman (2000), p. 37
    approx <- c(NA, NA, NA)
    if (n > 0){
        r <- round(n / 2 - (qnorm(1 - alpha / 2) * sqrt(n) / 2))
        s <- round(1 + n / 2 + (qnorm(1 - alpha / 2) * sqrt(n) / 2))
        approx <- c(m, v[r], v[s])
        }
    
    r <- data.frame(rbind(c(exact[2], exact[1], exact[3]), c(approx[2], approx[1], approx[3])))
    rownames(r) <- c("exact", "approximate")
    colnames(r) <- c("lower", "median", "upper")
    return(r)
}
