\name{confIntOddsRatio}
\alias{confIntOddsRatio}

\title{Compute confidence interval for an odds ratio}
\description{A methods to compute a confidence interval for an odds ratio
  is provided. The method is based on a Wald interval for the log odds
  ratio.   Used by \code{\link{confIntDiagnostic}}}
\usage{confIntOddsRatio(x, n, conf.level = 0.95)}
\arguments{
  \item{x}{vector of length 2, number of successes in each group.}
  \item{n}{vector of length 2, total number of trials in each group.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
  {A vector containing the odds ratio and the limits of the confidence interval.}
}

\references{DG Altman, D Machin, TN Bryant, MJ Gardner: Statistics with confidence, 2nd Edition, 2000, Chapter 7}

\author{Leonhard Held}
\seealso{\code{\link{confIntRiskRatio}}}

\examples{
x <- c(30, 50)
n <- c(100, 120)
confIntOddsRatio(x, n)
}

\keyword{htest}
