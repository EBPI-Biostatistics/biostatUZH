#' Sample size for McNemar test
#' 
#' Compute the sample size to test the null hypothesis of the McNemar test.
#' 
#' Given two samples of paired binary observations, e.g., results of some
#' positive/negative test on the same experimental units, we want to assess the
#' null hypothesis whether the probability of positives by the first method is
#' equal that of positives by the second method.
#' 
#' @param p1 Assumed marginal probability of the first row.
#' @param p2 Assumed marginal probability of the first column.
#' Has to be larger than p1 and 0.5.
#' @param sig.level Significance level.
#' @param power Power, i.e., 1 - probability of type II error.
#' @return Vector of min, max, and mid sample size as explained in
#' Lachenbruch (1982).
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Lachenbruch (1992). Sample Size for Studies based on McNemar's
#' test \emph{Statistics in Medicine}, \bold{11}, 1521--1525.
#' @keywords htest
#' @examples
#' 
#' # example from Lachenbruch (1982), Table II, first row
#' sampleSizeMcNemar(p1 = 0.8, p2 = 0.9, sig.level = 0.05, power = 0.9)
#' 
#' @export
sampleSizeMcNemar <- function(p1, p2, sig.level = 0.05, power = 0.9){

    stopifnot(is.numeric(p1), length(p1) == 1,
              is.finite(p1),
              0 <= p1, p1 <= 1,
              is.numeric(p2), length(p2) == 1,
              is.finite(p2),
              0 <= p2, p2 <= 1,
              p1 < p2, 
              0.5 < p2,
              is.numeric(sig.level), length(sig.level) == 1,
              is.finite(sig.level),
              0 < sig.level, sig.level < 1,
              is.numeric(power), length(power) == 1,
              is.finite(power),
              0 < power, power < 1)
        
    pp1 <- p1
    p1p <- p2
    p11 <- sort(seq(min(p1p, pp1), p1p + pp1 - 1, by = -10^-4))
    s <- (pp1 - p11) / (pp1 + p1p - 2 * p11)
    nl <- ceiling(0.25 * (qnorm(sig.level / 2) + qnorm(1 - power))^2 / (0.5 - s)^2  / (abs(pp1 + p1p - 2 * p11)))
    N <- nl[c(1, median(1:length(nl)), length(nl))]
    names(N) <- c("N_l min", "N_l mid", "N_l max")
    N
}
