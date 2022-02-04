\name{confIntRiskDiff}
\alias{confIntRiskDiff}

\title{Compute confidence interval for a risk difference}
\description{Two methods to compute a confidence interval for a risk difference
based on Wald and Wilson confidence intervals for the individual risks are provided.}
\usage{confIntRiskDiff(x, n, conf.level = 0.95)}
\arguments{
  \item{x}{vector of length 2, number of successes in each group.}
  \item{n}{vector of length 2, total number of trials in each group.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{rd}{Estimated risk difference.}
  \item{CIs}{Dataframe containing confidence intervals for the risk difference.}
}

\references{DG Altman, D Machin, TN Bryant, MJ Gardner: Statistics with confidence, 2nd Edition, 2000, Chapter 6}

\author{Leonhard Held}
\seealso{\code{\link{wilson}}}

\examples{
x <- c(30, 50)
n <- c(100, 120)
confIntRiskDiff(x, n)$CIs
}

\keyword{htest}
