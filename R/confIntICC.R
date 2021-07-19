#' Confidence intervals for intraclass correlation in interrater reliability
#' 
#' Compute, under suitable assumptions, confidence intervals for interrater
#' reliability for continuous measurements presented to two or more raters.
#' 
#' This function computes all the confidence intervals that are discussed in
#' Roussen et al. (2003). In applications, the interval under the "trained
#' rater" assumptions is often suitable.
#' 
#' @aliases confIntICC Aalpha lamAlpha rootFct
#' @param dat Data frame that contains the columns score, pat, rater.
#' @param conf.level Confidence level for confidence interval.
#' @param psi.re.0 2-d vector specifying the interval \eqn{[psi_0, psi_1]} on
#' p. 621 of Rousson et al. (2003).
#' @return A list containing: \item{icc(2, 1)}{ICC(2, 1): Intraclass
#' correlation from a two-random effects model.} \item{icc(3, 1)}{ICC(3, 1):
#' Intraclass correlation from a model with fixed rater effect.}
#' \item{psi_r/e}{The value \eqn{\psi_{r/e}} computed from the actual data.}
#' \item{ci.trained.rater}{Confidence interval under the trained rater
#' assumption, see Rousson et al. (2003), Section 4.}
#' \item{ci.low.asy.corr}{Lower bound of asymptotically exact confidence
#' interval, see Rousson et al. (2003), Section 3.}
#' \item{ci.low.fix.rater}{Lower bound of confidence interval under the fixed
#' rater assumption, see Rousson et al. (2003), Section 5.}
#' @note The function \code{\link{computeICCrater}} computes ICCs relying on a
#' mixed-model formulation, and is therefore able to handle unbalanced data. On
#' the contrary, the confidence intervals in the function
#' \code{\link{confIntICC}} are computed using sums of squares, and the data
#' must therefore be \emph{balanced}. See the example below.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @seealso \code{\link{computeICCrater}}
#' @references Rousson, V., Gasser, T., and Seifert, B. (2002). Assessing
#' intrarater, interrater and test-retest reliability of continuous
#' measurements. \emph{Statist. Med.}, \bold{21}, 3431--3446.
#' 
#' Rousson, V., Gasser, T., and Seifert, B. (2003). Confidence intervals for
#' intraclass correlation in inter-rater reliability. \emph{Scand. J.
#' Statist.}, \bold{30}, 617--624.
#' @keywords htest models
#' @examples
#' 
#' 
#' ## Generate dataset. Data must be balanced!
#' set.seed(1977)
#' n <- 40
#' r1 <- round(runif(n, 1, 20))
#' dat <- data.frame(
#'     "score" = c(r1, r1 + abs(round(rnorm(n, 1, 3))), 
#'         r1 + abs(round(rnorm(n, 1, 3)))), 
#'     "pat" = rep(c(1:n), 3),
#'     "rater" = rep(1:3, each = n)
#' )
#' confIntICC(dat, conf.level = 0.95, psi.re.0 = c(0, 1))
#' 
#' if (requireNamespace("lme4")) {
#'     computeICCrater(dat)
#' }
#' 
confIntICC <- function(dat, conf.level = 0.95, psi.re.0 = c(0, 1)){

## We compute the ICC(2, 1): Each subject is measured by each rater,
## and raters are considered representative of a larger population of
## similar raters. Reliability calculated from a single measurement.
##
## We also provide confidence intervals discussed in Rousson et al (2003),
## Scand J Statist 30, 617-624.
##
## Input:
##     - A data.frame "dat" containing the columns score, pat, rater.
##     - Confidence level alpha for confidence interval.
##     - psi.re.0: interval [psi_0, psi_1].
##
## Output:
##     - 
##
## Kaspar Rufibach, June 25, 2008
##

alpha <- 1 - conf.level

## sort by rater, pat
dat <- dat[order(dat$rater, dat$pat), ]

## compute ICC using sum of squares
## see Rousson et al (2002), Stat in Med
rj <- unique(dat$rater)
n <- sum(dat$rater == rj[1])
d <- length(rj)
Yij <- dat$score
Yi <- as.vector(unlist(lapply(split(Yij, dat$pat), mean)))
Yj <- as.vector(unlist(lapply(split(Yij, dat$rater), mean)))
Y <- mean(Yij)

MSs <- d / (n - 1) * sum((Yi - Y) ^ 2)
MSr <- n / (d - 1) * sum((Yj - Y) ^ 2)

tmp <- 0
for (j in 1:d){tmp <- tmp + ((Yij - Yi)[dat$rater == rj[j]] - Yj[j] + Y) ^ 2}
MSe <- 1 / ((d - 1) * (n - 1)) * sum(tmp)
sig2.ms <- c((MSs - MSe) / d, (MSr - MSe) / n, MSe)

icc21 <- sig2.ms[1] / sum(sig2.ms)
icc31 <- sig2.ms[1] / sum(sig2.ms[c(1, 3)])

## compute confidence intervals
# define grid of psi_{s/e}'s
psi.se <- sig2.ms[1] / sig2.ms[3]
psi.re <- sig2.ms[2] / sig2.ms[3]

# asymptotically correct CI, p. 620
Ln <- (MSs - MSe) / (MSs + d * (d - 1) * MSr / (n * qchisq(alpha / 2, df = d - 1)) + (d - 1) * MSe)

# population of trained raters: L_n on p. 621
Ltr <- Aalpha(alpha / 2, n, d, MSs, MSe) / (Aalpha(alpha / 2, n, d, MSs, MSe) + psi.re.0[2] + 1)
Utr <- Aalpha(1 - alpha / 2, n, d, MSs, MSe) / (Aalpha(1 - alpha / 2, n, d, MSs, MSe) + psi.re.0[1] + 1)

# model with fixed rater effect: interval on p. 622
B1 <- lamAlpha(alpha / 2, n, d, MSs, MSe) / ((d - 1) * n)
B2 <- lamAlpha(1 - alpha / 2, n, d, MSs, MSe) / ((d - 1) * n)

# lower bound for rho.til at level (1 - alpha)
L.til <- Aalpha(alpha / 4, n, d, MSs, MSe) / (Aalpha(alpha / 4, n, d, MSs, MSe) + B2 + 1)

# generate output
res <- list("ICC(2, 1)" = icc21, "ICC(3, 1)" = icc31, "psi_r/e" = psi.re, "ci.trained.rater" = c(Ltr, Utr), "ci.low.asy.corr" = Ln, "ci.low.fix.rater" = L.til)
return(res)
}
