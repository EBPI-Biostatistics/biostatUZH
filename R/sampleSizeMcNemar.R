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
#' # example from Lachenbruch (1992), Table II, first row
#' sampleSizeMcNemar(p1 = 0.8, p2 = 0.9, sig.level = 0.05, power = 0.9)
#'
#' @export
sampleSizeMcNemar <- function(
    p1, p2, # must be given for Lachenbruch method
    p01, p10, # must be given for Connor method
    sig.level = 0.05,
    power = 0.9,
    method = c("lb", "co"),
    cc = TRUE
) {

    # Match method argument
    method <- match.arg(method)

    # Input checks
    if (method == "lb") {
        if (missing(p1) || missing(p2)) {
            stop("When specifying 'method = \"lb\"', arguments p1 and p2 must be specified.")
        }
        if (!missing(p01) || !missing(p10)) {
            message("Arguments 'p01' and 'p10' are ignored when 'method = \"lb\"'.")
        }

        stopifnot(
            is.numeric(p1), length(p1) == 1,
            is.finite(p1),
            0 <= p1, p1 <= 1,
            is.numeric(p2), length(p2) == 1,
            is.finite(p2),
            0 <= p2, p2 <= 1,
            p1 < p2,
            0.5 < p2
        )

        probs <- c(p1, p2)

    } else { #if (method == "co") {
        if (missing(p01) || missing(p10))
            stop("When specifying 'method = \"co\"', arguments p01 and p10 must be specified.")
        if (!missing(p1) || !missing(p2))
            message("Arguments 'p1' and 'p2' are ignored when 'method = \"co\"'.")

        stopifnot(
            is.numeric(p01), length(p01) == 1,
            is.finite(p01),
            0 <= p01, p01 <= 1,
            is.numeric(p10), length(p10) == 1,
            is.finite(p10),
            0 <= p10, p10 <= 1,
        )

        probs <- c(p01, p10)

    }

    stopifnot(
        is.numeric(sig.level), length(sig.level) == 1,
        is.finite(sig.level),
        0 < sig.level, sig.level < 1,
        is.numeric(power), length(power) == 1,
        is.finite(power),
        0 < power, power < 1
    )

    # put arguments in a list
    args <- list(
        probs[1],
        probs[2],
        sig.level = sig.level,
        power = power
    )

    # Determine the desired function
    f <- switch(
        method,
        "lb" = "sampleSizeLachenbruch",
        "co" = "sampleSizeConnor"
    )

    ceiling(do.call(f, args))
}

#' Comupute the sample size based on Lachenbruch's method
#'
#' @description This function implements the sample size calculation
#' according to Lachenbruch(1992).
#' @noRd
sampleSizeLachenbruch <- function(
    p1,
    p2,
    sig.level,
    power
) {

    pp1 <- p1
    p1p <- p2
    p11 <- sort(seq(min(p1p, pp1), p1p + pp1 - 1, by = -10^-4))
    s <- (pp1 - p11) / (pp1 + p1p - 2 * p11)
    nl <- 0.25 * (qnorm(sig.level / 2) + qnorm(1 - power))^2 /
        (0.5 - s)^2  / (abs(pp1 + p1p - 2 * p11))
    N <- nl[c(1, median(seq_along(nl)), length(nl))]
    names(N) <- c("N_l min", "N_l mid", "N_l max")
    N
}

#' Comupute the sample size based on Connor's method
#'
#' @description This function implements the sample size calculation
#' according to Connor(1987).
#' @noRd
sampleSizeConnor <- function(
    p01, p10,
    sig.level,
    power,
    cc
) {

    # psi refers to the probability of discordance
    psi <- p01 + p10
    # This is the difference between probabilities of discordance
    delta <- p01 - p10
    n <- (
        (
            qnorm(sig.level / 2) * sqrt(psi) +
                qnorm(1 - power) * sqrt(psi - delta^2)
        ) / delta
    )^2
    # Continuity correction
    if (cc) {
        n <- n + 1 / abs(delta)
    }
    n
}

sampleSizeConnor(
    p01 = 0.25, p10 = 0.15, sig.level = 0.05, power = 0.8, cc = FALSE
)
