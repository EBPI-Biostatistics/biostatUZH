#' Confidence interval for binomial proportions
#' 
#' Compute a confidence interval for binomial proportions using several
#' asymptotic and exact methods.
#' 
#' 
#' @aliases confIntProportion wald wilson agresti jeffreys clopperPearson
#' @param x Number of successes.
#' @param n Total number of trials.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return \code{confIntProportion} returns a data.frame with confidence
#' intervals from the Wald, Wilson, Agresti, Jeffreys, and Clopper-Pearson
#' methods.
#' 
#' \code{wald} returns the Wald confidence interval.
#' 
#' \code{wilson} returns the Wilson confidence interval.
#' 
#' \code{agresti} returns the Agresti confidence interval.
#' 
#' \code{jeffreys} returns the Jeffreys confidence interval.
#' 
#' \code{clopperPearson} returns the Clopper-Pearson confidence interval.
#' @author Kaspar Rufibach, Leonhard Held
#' @seealso Functions for some of the intervals provided here are available in
#' \pkg{Hmisc}; see the examples.
#' @references All the intervals provided in these functions are compared in:
#' 
#' Brown, L.D., Cai, T.T., DasGupta, A. (2001). Interval Estimation for a
#' Binomial Proportion. \emph{Statistical Science}, \bold{16(2)}, 101--133.
#' @examples
#' 
#' ## Calculate confidence bounds for a binomial parameter by different methods.
#' x <- 50
#' n <- 100
#' ci <- confIntProportion(x, n)$CIs
#' ci
#' 
#' plot(0, 0, type = 'n', ylim = c(0, 7), xlim = c(0, 1), xlab = 'p',
#'      ylab = '', yaxt = 'n')
#' for(i in 1:5)
#'     lines(ci[i, 2:3], c(i, i))
#' text(0.5, 0.85, 'wald')
#' text(0.5, 1.85, 'wilson')
#' text(0.5, 2.85, 'agresti')
#' text(0.5, 3.85, 'jeffreys')
#' text(0.5, 4.85, 'clopper')
#' 
#' ## compare intervals to those received by the function binconf in Hmisc:
#' if (require("Hmisc")) {
#'     binconf(x, n, method = "asymptotic")   # Wald
#'     binconf(x, n, method = "wilson")       # Wilson
#'     binconf(x, n, method = "exact")        # Clopper-Pearson
#' }
#' 
confIntProportion <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1)

    res <- data.frame(matrix(NA, ncol = 3))
    colnames(res) <- c("type", "lower", "upper")

    res[1, 2:3] <- wald(x, n, conf.level = conf.level)[c(1, 3)]
    res[2, 2:3] <- wilson(x, n, conf.level = conf.level)[c(1, 3)]
    res[3, 2:3] <- agresti(x, n, conf.level = conf.level)[c(1, 3)]
    res[4, 2:3] <- jeffreys(x, n, conf.level = conf.level)[c(1, 3)]
    res[5, 2:3] <- clopperPearson(x, n, conf.level = conf.level)[c(1, 3)]

    res[, 1] <- c("Wald", "Wilson", "Agresti", "Jeffreys", "ClopperPearson")

    res <- list("p" = x / n, "CIs" = res)
    return(res)
}


#' @rdname confIntProportion
#' @return \code{wald} returns the Wald confidence interval.
#' @export 
wald <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1)
    q <- qnorm(p=(1+conf.level)/2)
    pi <- x/n         
    limits <- pi + c(-1, 1)*q*sqrt(pi*(1-pi)/n)
    res <- c("lower" = max(0, limits[1]), "prop" = pi, "upper" = min(1, limits[2]))
    return(res)
}



#' @rdname confIntProportion
#' @return \code{wilson} returns the Wilson confidence interval.
#' @export 
wilson <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1) # n=0 yields NaN results, but allow for summaryROC()
    q <- qnorm(p=(1+conf.level)/2)
    q2 <- q^2
    prop <- x/n
    mid <- (x+q2/2)/(n+q2)
    factor <- (q*sqrt(n))/(n+q2)*sqrt(prop*(1-prop)+q2/(4*n))
    limits <- mid + c(-1,1)*factor
    res <- c("lower"=limits[1], "prop"=prop, "upper"=limits[2])
    return(res)
}

#' @rdname confIntProportion
#' @return \code{agresti} returns the Agresti confidence interval.
#' @export 
agresti <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1)
    k <- qnorm(p = (conf.level + 1) / 2)
    ptilde <- (x + 2) / (n + 4)
    z <- abs(k)
    stderr <- sqrt(ptilde * (1 - ptilde) / (n + 4))
    ll <- max(ptilde - z * stderr, 0)
    ul <- min(ptilde + z * stderr, 1)
    res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
    
    return(res)
}

#' @rdname confIntProportion
#' @return \code{jeffreys} returns the Jeffreys confidence interval.
#' @export 
jeffreys <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1)
    q <- (1-conf.level)/2
    alpha <- x + 0.5
    beta <- n - x + 0.5
    pihat <- qbeta(0.5, alpha, beta)
    limits <- qbeta(c(q, 1-q), alpha, beta)
    res <- c("lower" = limits[1], "pihat" = pihat, "upper" = limits[2])
    return(res)
}

#' @rdname confIntProportion
#' @return \code{clopperPearson} returns the Clopper-Pearson confidence interval.
#' @export 
clopperPearson <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), is.wholenumber(n), x<=n, n>=1,
              0 < conf.level, conf.level<1)
    a <- 1 - conf.level
    if (x == 0){
        ll <- 0
        ul <- 1 - (a / 2) ^ (1 / n)}
    else if (x == n){
        ll <- (a / 2) ^ (1 / n)
        ul <- 1
    }
    else {
        ll <- 1/(1 + (n - x + 1) / (x * qf(a / 2, 2 * x, 2 * (n - x + 1))))
        ul <- 1/(1 + (n - x) / ((x + 1) * qf(1 - a / 2, 2 * (x + 1), 2 * (n - x))))
    }
    res <- c("lower" = ll, "prop" = x / n, "upper" = ul)
    return(res)
} 
