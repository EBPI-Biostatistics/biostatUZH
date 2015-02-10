\name{computeICCrater}
\alias{computeICCrater}
\title{Compute intraclass correlation coefficients ICC(2, 1) and ICC(3, 1)}
\description{Interrater reliability for continuous measurements presented to two or more raters.}
\usage{computeICCrater(dat)}
\arguments{\item{dat}{Data frame that contains the columns score, pat, rater.}}
\details{Interrater reliability, or more specifically an intraclass correlation coefficient (ICC), is computed when 
the same continuous measurements are presented to two or more raters. Depending on the underlying model, several 
different ICCs can be computed. Here, we provide ICC(2, 1) which is derived from a regression model with random
effects for rater and patients, and ICC(3, 1) where raters are considered to be fixed.}
\value{A list containing the following elements:
\item{sig.pat}{Variance of random subject effect.}
\item{sig.rater}{Variance of the random rater effect.}
\item{sig.res}{Residual variance.}
\item{icc(2, 1)}{ICC(2, 1): Intraclass correlation from a two-random effects model.}
\item{icc(3, 1)}{ICC(3, 1): Intraclass correlation from a model with fixed rater effect.}
}

\note{Since this function relies on a mixed-model formulation using
  package \pkg{lme4} to compute ICCs, it can also handle \emph{unbalanced}
  data.}


\references{
Shrout, P.E. and Fleiss, J.L. (1979).
Intraclass correlations: uses in assessing rater reliability.
\emph{Psychological Bulletin}, \bold{36}, 420--428.

Rousson, V., Gasser, T., and Seifert, B. (2002).
Assessing intrarater, interrater and test-retest reliability of continuous measurements.
\emph{Statist. Med.}, \bold{21}, 3431--3446.
}
\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}
\seealso{\code{\link{confIntICC}}}

\examples{
if (requireNamespace("lme4")) {
  ## generate dataset (balanced data)
  set.seed(1977)
  n <- 40
  r1 <- round(runif(n, 1, 20))
  dat <- data.frame(
      "score" = c(r1, r1 + abs(round(rnorm(n, 1, 3))),
          r1 + abs(round(rnorm(n, 1, 3)))), 
      "pat" = rep(c(1:n), 3),
      "rater" = rep(1:3, each = n))
  computeICCrater(dat)

  ## also works for unbalanced data
  dat2 <- dat[sort(sample(1:(3 * n))[1:100]), ]
  computeICCrater(dat2)
}
}

\keyword{htest}
\keyword{models}

