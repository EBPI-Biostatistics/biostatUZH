\name{confIntICC}
\alias{confIntICC}
\alias{Aalpha}
\alias{lamAlpha}
\alias{rootFct}
\title{Confidence intervals for intraclass correlation in interrater reliability}
\description{Compute, under suitable assumptions, confidence intervals for interrater reliability for continuous 
measurements presented to two or more raters.}
\usage{
confIntICC(dat, conf.level = 0.95, psi.re.0 = c(0, 1))
}
\arguments{
    \item{dat}{Data frame that contains the columns score, pat, rater.}  
    \item{conf.level}{Confidence level for confidence interval.}
    \item{psi.re.0}{2-d vector specifying the interval \eqn{[psi_0, psi_1]} on p. 621 of Rousson et al. (2003).}
}
\details{This function computes all the confidence intervals that are discussed in Roussen et al. (2003).
In applications, the interval under the "trained rater" assumptions is often suitable.}
\value{A list containing:
\item{icc(2, 1)}{ICC(2, 1): Intraclass correlation from a two-random effects model.}
\item{icc(3, 1)}{ICC(3, 1): Intraclass correlation from a model with fixed rater effect.}
\item{psi_r/e}{The value \eqn{\psi_{r/e}} computed from the actual data.}
\item{ci.trained.rater}{Confidence interval under the trained rater assumption, see Rousson et al. (2003), Section 4.}
\item{ci.low.asy.corr}{Lower bound of asymptotically exact confidence interval, see Rousson et al. (2003), Section 3.}
\item{ci.low.fix.rater}{Lower bound of confidence interval under the fixed rater assumption, see Rousson et al. (2003), Section 5.}
}

\note{The function \code{\link{computeICCrater}} computes ICCs relying on a mixed-model formulation, and is therefore
able to handle unbalanced data. On the contrary, the confidence intervals in the function \code{\link{confIntICC}} are
computed using sums of squares, and the data must therefore be \emph{balanced}. See the example below.}

\references{
Rousson, V., Gasser, T., and Seifert, B. (2002).
Assessing intrarater, interrater and test-retest reliability of continuous measurements.
\emph{Statist. Med.}, \bold{21}, 3431--3446.

Rousson, V., Gasser, T., and Seifert, B. (2003).
Confidence intervals for intraclass correlation in inter-rater reliability.
\emph{Scand. J. Statist.}, \bold{30}, 617--624.
}
\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com} }
\seealso{\code{\link{computeICCrater}}}
\examples{

## Generate dataset. Data must be balanced!
set.seed(1977)
n <- 40
r1 <- round(runif(n, 1, 20))
dat <- data.frame(
    "score" = c(r1, r1 + abs(round(rnorm(n, 1, 3))), 
        r1 + abs(round(rnorm(n, 1, 3)))), 
    "pat" = rep(c(1:n), 3),
    "rater" = rep(1:3, each = n)
)
confIntICC(dat, conf.level = 0.95, psi.re.0 = c(0, 1))

if (requireNamespace("lme4")) {
    computeICCrater(dat)
}
}

\keyword{htest}
\keyword{models}
