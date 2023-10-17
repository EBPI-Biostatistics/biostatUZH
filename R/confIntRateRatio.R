#' @name confIntRateRatio
#'
#' @title Compute confidence interval for a rate ratio
#'
#' @description A method to compute a confidence interval for a rate ratio
#'     is provided. The method is based on a Wald interval for the log rate
#'     ratio.
#'
#' @param x A vector of length 2 containing the number of events in each group.
#' @param t A numeric vector of length 2 containing the total observation time
#'     in each group.
#' @param conf.level A numeric vector of length 1 containing the confidence
#'     level for the confidence interval.
#'
#' @examples
#'     x <- c(30, 50)
#'     t <- c(100, 120)
#'     confIntRateRatio(x, t)
#'
#' @return A named numeric vector of length 3 containing the rate ratio and
#'     the limits of the confidence interval.
#'
#' @references Held, L., Rufibach, K. and Seifert, B. (2013)
#'     \emph{Medizinische Statistik - Konzepte, Methoden, Anwendungen}.
#'     Section 8.3.
#'
#' @author Leonhard Held
#'
#' @seealso \code{link[biostatUZH]{confIntRateDiff}}
#'
#' @keywords htest
#'
#' @export
confIntRateRatio <- function(x, t, conf.level = 0.95) {

    stopifnot(
        is.numeric(x),
        is.numeric(t),
        is.finite(x),
        is.finite(t),
        length(x) == 2,
        length(t) == 2,
        is.wholenumber(x),
        all(t > 0),
        all(x > 0),
        is.numeric(conf.level),
        is.finite(conf.level),
        length(conf.level) == 1,
        conf.level < 1,
        conf.level > 0
    )

    Rate <- x / t
    RateRatio <- Rate[1] / Rate[2]
    se.log.RateRatio <- sqrt(sum(1 / x))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * se.log.RateRatio)
    wald.lower <- RateRatio / EF
    wald.upper <- RateRatio * EF

    c(
        "lower" = wald.lower,
        "Rate Ratio" = RateRatio,
        "upper" = wald.upper
    )
}
