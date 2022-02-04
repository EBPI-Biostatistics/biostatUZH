\name{mantelHaenszel}
\alias{mantelHaenszel}
\title{Cochran-Mantel-Haenszel Chi-Squared Test for Count Data}
\description{Compute Cochran-Mantel-Haenszel chi-squared test of the null that two nominal variables are conditionally 
independent in each stratum, assuming that there is no three-way interaction.}
\usage{mantelHaenszel(exposure, outcome, stratum)}
\arguments{
  \item{exposure}{Binary variable coding whether patient was exposed (1) or not (0).}
  \item{outcome}{Binary variable coding outcome of patient.}
  \item{stratum}{Nominal variable containing information about matching, i.e. strata.}
}

\value{
\item{tab}{Table that counts numbers of strata for each case-control combination.}
\item{test.stat}{Test statistic for \eqn{chi^2} test.}
\item{p.val}{\eqn{p}-value of \eqn{chi^2} test.}
}

\references{Agresti, A. (2002). \emph{Categorical data analysis}. New York: Wiley.}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com} }

\seealso{Similar functionality is provided in \code{\link{mantelhaen.test}} and \code{clogit} in \pkg{survival}.
See the examples below for a comparison.}

\examples{

# generate data
set.seed(1977)
exposure <- rep(c(1, 0, 0, 0, 0), 41)
outcome <- sample(c(rep(1, 62), rep(0, 5 * 41 - 62)))
strata <- rep(1:41, each = 5)

# via conditional logistic regression
logreg <- clogit(outcome ~ exposure + strata(strata), method = "approximate")
summary(logreg)

# R function in library 'stats'
mh <- mantelhaen.test(x = outcome, y = exposure, z = strata)

# this function
mH <- mantelHaenszel(exposure, outcome, strata)

# compare p-values
summary(logreg)$coef[5]
mh$p.value
mH$p.val
}

\keyword{htest}
\keyword{survival}
