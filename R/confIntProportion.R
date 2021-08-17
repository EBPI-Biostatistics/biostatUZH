#' Confidence interval for binomial proportions
#' 
#' Compute a confidence interval for binomial proportions using several
#' asymptotic and exact methods.
#' 
#' @aliases confIntProportion wald wilson agresti jeffreys clopperPearson
#' @param x Number of successes.
#' @param n Total number of trials.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return \code{confIntProportion} returns a data.frame with confidence
#' intervals from the Wald, Wilson, Agresti, Jeffreys, and Clopper-Pearson
#' methods.
#' @author Kaspar Rufibach, Leonhard Held, Florian Gerber
#' @seealso Functions for some of the intervals provided here are available in
#' \pkg{Hmisc}. \code{\link{confIntIndependentProportion}}
#' @references All the intervals provided in these functions are compared in:
#' 
#' Brown, L.D., Cai, T.T., DasGupta, A. (2001). Interval Estimation for a
#' Binomial Proportion. \emph{Statistical Science}, \bold{16(2)}, 101--133.
#' @examples
#' ## Calculate confidence bounds for a binomial parameter using different methods.
#' x <- 50
#' n <- 100
#' ci <- confIntProportion(x = x, n = n)$CIs
#' ci
#' 
#' plot(0, 0, type = "n", ylim = c(0, 7), xlim = c(0, 1), xlab = "p",
#'      ylab = "", yaxt = "n")
#' for(i in 1:5)
#'     lines(ci[i, 2:3], c(i, i))
#' text(0.5, 0.85, "wald")
#' text(0.5, 1.85, "wilson")
#' text(0.5, 2.85, "agresti")
#' text(0.5, 3.85, "jeffreys")
#' text(0.5, 4.85, "clopper")
#' 
#' ## compare intervals to those received by the function binconf in Hmisc:
#' if (require("Hmisc")) {
#'     binconf(x, n, method = "asymptotic")   # Wald
#'     binconf(x, n, method = "wilson")       # Wilson
#'     binconf(x, n, method = "exact")        # Clopper-Pearson
#' }
#' 
#' @export
confIntProportion <- function(x, n, conf.level = 0.95)
{

    ## input tests are done in wald()
    
    res <- data.frame(matrix(NA, ncol = 3))
    colnames(res) <- c("type", "lower", "upper")

    res[1, 2:3] <- wald(x, n, conf.level = conf.level)[c(1, 3)]
    res[2, 2:3] <- wilson(x, n, conf.level = conf.level)[c(1, 3)]
    res[3, 2:3] <- agresti(x, n, conf.level = conf.level)[c(1, 3)]
    res[4, 2:3] <- jeffreys(x, n, conf.level = conf.level)[c(1, 3)]
    res[5, 2:3] <- clopperPearson(x, n, conf.level = conf.level)[c(1, 3)]

    res[, 1] <- c("Wald", "Wilson", "Agresti", "Jeffreys", "ClopperPearson")

    list("p" = x / n, "CIs" = res)
}


#' @rdname confIntProportion
#' @return \code{wald} returns the Wald confidence interval.
#' @export 
wald <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.numeric(x), length(x) == 1,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 1,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 1,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)

    q <- qnorm(p = (1 + conf.level) / 2)
    prop <- x / n         
    limits <- prop + c(-1, 1) * q * sqrt(prop * (1 - prop) / n)
    c("lower" = max(0, limits[1]), "prop" = prop, "upper" = min(1, limits[2]))

}



#' @rdname confIntProportion
#' @return \code{wilson} returns the Wilson confidence interval.
#' @export 
wilson <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.numeric(x), length(x) == 1,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 1,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 0, # n=0 yields NaN results, but allow for summaryROC()
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    q <- qnorm(p = (1 + conf.level) / 2)
    q2 <- q^2
    prop <- x / n
    mid <- (x + q2 / 2) / (n + q2)
    factor <- (q * sqrt(n)) / (n + q2) * sqrt(prop * (1 - prop) + q2 / (4 * n))
    limits <- mid + c(-1, 1) * factor
    c("lower" = limits[1], "prop" = prop, "upper" = limits[2])
}

#' @rdname confIntProportion
#' @return \code{agresti} returns the Agresti confidence interval.
#' @export 
agresti <- function(x, n, conf.level = 0.95)
{
        stopifnot(is.numeric(x), length(x) == 1,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 1,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 1,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    k <- qnorm(p = (conf.level + 1) / 2)
    ptilde <- (x + 2) / (n + 4)
    z <- abs(k)
    stderr <- sqrt(ptilde * (1 - ptilde) / (n + 4))
    ll <- max(ptilde - z * stderr, 0)
    ul <- min(ptilde + z * stderr, 1)
    c("lower" = ll, "prop" = x / n, "upper" = ul)
}

#' @rdname confIntProportion
#' @return \code{jeffreys} returns the Jeffreys confidence interval.
#' @export 
jeffreys <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.numeric(x), length(x) == 1,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 1,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 1,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    q <- (1 - conf.level) / 2
    alpha <- x + 0.5
    beta <- n - x + 0.5
    pihat <- qbeta(0.5, alpha, beta)
    limits <- qbeta(c(q, 1 - q), alpha, beta)
    c("lower" = limits[1], "pihat" = pihat, "upper" = limits[2])
}

#' @rdname confIntProportion
#' @return \code{clopperPearson} returns the Clopper-Pearson confidence interval.
#' @export 
clopperPearson <- function(x, n, conf.level = 0.95)
{
    stopifnot(is.numeric(x), length(x) == 1,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 1,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 1,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    a <- 1 - conf.level
    if (x == 0){
        ll <- 0
        ul <- 1 - (a / 2) ^ (1 / n)
    } else if (x == n){
        ll <- (a / 2) ^ (1 / n)
        ul <- 1
    } else {
        ll <- 1/(1 + (n - x + 1) / (x * qf(a / 2, 2 * x, 2 * (n - x + 1))))
        ul <- 1/(1 + (n - x) / ((x + 1) * qf(1 - a / 2, 2 * (x + 1), 2 * (n - x))))
    }
    c("lower" = ll, "prop" = x / n, "upper" = ul)
} 


#' Confidence interval for the risk difference of two independent samples
#' 
#' Compute confidence interval for the risk difference of two independent
#' samples based on individual Wilson intervals using Newcombe's method.
#' 
#' 
#' @aliases confIntIndependentProportion Newcombe
#' @param x Vector with two entries, the successes in the two groups.
#' @param n Vector with two entries, the number of trials.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries:
#' \item{p1}{Estimated proportion in first sample.}
#' \item{p2}{Estimated proportion in second sample.}
#' \item{d}{Estimated difference \eqn{p_1 - p_2} of proportions.}
#' \item{newcombeCI}{Confidence interval for the difference of independent
#' proportions, computed according to Newcombe's method.}
#' \item{waldCI}{Wald confidence interval for the difference of independent
#' proportions.}
#' @author Leonhard Held
#' @seealso \code{\link{confIntPairedProportion}}, \code{\link{confIntProportion}},
#' \code{\link[stats]{fisher.test}}, \code{\link[exact2x2]{exact2x2}}
#' @references The Newcombe interval is introduced in
#' 
#' Newcombe, R.G. (1998). Interval estimation for the difference between
#' independent proportions: Comparison of eleven methods. \emph{Stat. Med.},
#' \bold{17}, 873--890.
#' 
#' A worked out example can be found in
#' 
#' Altman, D.G., Machin, D., Bryant, T.N., Gardner, M.J. (2000). Statistics
#' with confidence (p. 49). University Press Belfast.
#' @keywords univar htest
#' @examples
#' 
#' # Example from Significance (2010), 7(4), p. 146, "Untimely ripped ?"
#' n <- c(1515, 108000)
#' x <- c(0, 100)
#' confIntIndependentProportion(x = x, n = n)
#' 
#' # Fisher p-value
#' tab <- matrix(data = c(x, n - x), nrow = 2, ncol = 2)
#' fisher <- fisher.test(x = tab)
#' fisher$p.value
#' 
#' # exact tests
#' if(require("exact2x2")){
#'     exact2x2(x = tab, tsmethod = "minlike")$p.value
#'     exact2x2(x = tab, tsmethod = "central")$p.value
#'     exact2x2(x = tab, tsmethod = "blaker")$p.value
#' }
#' 
#' @export
confIntIndependentProportion <- function(x, n, conf.level = 0.95)
{
    ## Calculates confidence interval for risk difference of two independent
    ## samples based on individual Wilson intervals
    stopifnot(is.numeric(x), length(x) == 2,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 2,
              is.finite(n), is.wholenumber(n),
              x <= n, n >= 1,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    alpha <- 1 - conf.level

    p1 <- x[1] / n[1]
    p2 <- x[2] / n[2]
    D <- p1 - p2

    ## Newcombe interval
    w1 <- wilson(x = x[1], n = n[1], conf.level = conf.level)
    w2 <- wilson(x = x[2], n = n[2], conf.level = conf.level)
    D.lower <- D - sqrt((p1 - w1[1]) ^ 2 + (p2 - w2[3]) ^ 2)
    D.upper <- D + sqrt((p1 - w1[3]) ^ 2 + (p2 - w2[1]) ^ 2)
    newcombe <- c(D.lower, D.upper)

    ## Wald interval
    se.D <- sqrt(p1 * (1 - p1) / n[1] + p2 * (1 - p2) / n[2])
    factor <- qnorm(1 - alpha / 2)
    D.lower <- D - factor * se.D
    D.upper <- D + factor * se.D
    wald <- c(D.lower, D.upper)

    list("p1" = p1, "p2" = p2, "d" = as.numeric(D),
         "newcombeCI" = as.numeric(newcombe),
         "waldCI" = as.numeric(wald))
}


#' Confidence interval for the difference of paired binomial
#' proportions using Newcombe's method
#' 
#' Compute confidence interval for the difference of paired binomial
#' proportions using Newcombe's method.
#' 
#' 
#' @aliases confIntPairedProportion
#' @param x A two-dimensional contingency table in matrix form.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries:
#' \item{p1}{Estimated proportion \eqn{p_{1+}}.}
#' \item{p2}{Estimated proportion \eqn{p_{+1}}.}
#' \item{newcombeCI}{Confidence interval for the difference of paired
#' proportions, computed according to Newcombe (1998).}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references The Newcombe interval is introduced in
#' 
#' Newcombe, R.G. (1998). Improved confidence intervals for the difference
#' between binomial proportions based on paired data. \emph{Stat. Med.},
#' \bold{17}, 2635--2650.
#' 
#' A worked out example can be found in
#' 
#' Altman, D.G., Machin, D., Bryant, T.N., Gardner, M.J. (2000). Statistics
#' with confidence. University Press Belfast.
#' @keywords univar htest
#' @seealso \code{\link{confIntIndependentProportion}}, \code{\link{wilson}},
#' \code{\link[exact2x2]{mcnemar.exact}}
#' @examples
#' 
#' # Calculate confidence interval for the example in Altman et al (2000), Table 6.2
#' altman62 <- rbind(c(14, 5), c(0, 22))
#' confIntPairedProportion(x = altman62)
#' 
#' # exact McNemar test
#' if(require("exact2x2")){
#'     mcnemar.exact(altman62)
#'     mcnemar.test(altman62)
#' }
#' 
#' @export
confIntPairedProportion <- function(x, conf.level = 0.95){
   
    ## confidence interval for a paired proportion, according to Newcombe (1998)
    ## see Altman (2000), "Statistics with confidence", p. 52 for a worked out example

    stopifnot(is.numeric(x), is.matrix(x), dim(x) == c(2, 2),
              is.finite(x), is.wholenumber(x),
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    r <- x[1, 1]
    s <- x[1, 2]
    t <- x[2, 1]
    u <- x[2, 2]
    n <- r + s + t + u
    
    p1 <- (r + s) / n
    p2 <- (r + t) / n
    D <- p1 - p2
    
    ## compute Wilson CIs for p1 and p2
    ci1 <- wilson(x = r + s, n = n, conf.level = conf.level)
    ci2 <- wilson(x = r + t, n = n, conf.level = conf.level)
    l1 <- ci1[1]
    u1 <- ci1[3]
    l2 <- ci2[1]
    u2 <- ci2[3]

    phi <- 0
    if (max(r + s, t + u, r + t, s + u) > 0){
        A <- (r + s) * (t + u) * (r + t) * (s + u)
        B <- r * u - s * t
        C <- 0
        if (B > n / 2){
            C <- B - n / 2
        }
        if (B < 0){
            C <- B
        }
        phi <- C / sqrt(A)    
    }
    
    newcombe1 <- D - sqrt((p1 - l1) ^ 2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2) ^ 2)
    newcombe2 <- D + sqrt((p2 - l2) ^ 2 - 2 * phi * (p2 - l2) * (u1 - p1) + (u1 - p1) ^ 2)
    
    newcombe <- c(newcombe1, newcombe2)
    list("p1" = p1, "p2" = p2, "diff" = D, "newcombeCI" = as.numeric(newcombe))
}

