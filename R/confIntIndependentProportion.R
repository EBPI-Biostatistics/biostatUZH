#' Compute confidence interval for the risk difference of two independent
#' samples
#' 
#' Compute confidence interval for the risk difference of two independent
#' samples based on individual Wilson intervals using Newcombe's method.
#' 
#' 
#' @aliases confIntIndependentProportion Newcombe
#' @param x Vector with two entries, the successes in the two groups.
#' @param n Vector with two entries, the number of trials.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries: \item{p1}{Estimated proportion in first
#' sample.} \item{p2}{Estimated proportion in second sample.}
#' \item{d}{Estimated difference \eqn{p_1 - p_2} of proportions.}
#' \item{newcombeCI}{Confidence interval for the difference of independent
#' proportions, computed according to Newcombe's method.} \item{waldCI}{Wald
#' confidence interval for the difference of independent proportions.}
#' @author Leonhard Held
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
#' confIntIndependentProportion(x, n)
#' 
#' # Fisher p-values
#' t <- matrix(c(x, n - x), nrow = 2, ncol = 2)
#' f.t <- fisher.test(t)
#' f.t$p.value
#' 
#' \dontrun{
#' # exact test
#' library("exact2x2")
#' exact2x2(t, tsmethod = "minlike")$p.value
#' exact2x2(t, tsmethod = "central")$p.value
#' exact2x2(t, tsmethod = "blaker")$p.value
#' }
#' 
#' @export
confIntIndependentProportion <- function(x, n, conf.level = 0.95)
{
    ## Calculates confidence interval for risk difference of two independent
    ## samples based on individual Wilson intervals

    alpha <- 1 - conf.level

    p1 <- x[1] / n[1]
    p2 <- x[2] / n[2]
    D <- p1 - p2

    ## Newcombe interval
    w1 <- wilson(x[1], n[1], conf.level=conf.level)
    w2 <- wilson(x[2], n[2], conf.level=conf.level)
    D.lower <- D - sqrt((p1 - w1[1]) ^ 2 + (p2 - w2[3]) ^ 2)
    D.upper <- D + sqrt((p1 - w1[3]) ^ 2 + (p2 - w2[1]) ^ 2)
    newcombe <- c(D.lower, D.upper)

    ## Wald interval
    se.D <- sqrt(p1 * (1 - p1) / n[1] + p2 * (1 - p2) / n[2])
    factor <- qnorm(1 - alpha / 2)
    D.lower <- D - factor * se.D
    D.upper <- D + factor * se.D
    wald <- c(D.lower, D.upper)

    res <- list("p1" = p1, "p2" = p2, "d" = as.numeric(D),
                "newcombeCI" = as.numeric(newcombe),
                "waldCI" = as.numeric(wald))
    return(res)
}
