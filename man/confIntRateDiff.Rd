\name{confIntRateDiff}
\alias{confIntRateDiff}

\title{Compute confidence interval for a rate difference}
\description{Two methods to compute a confidence interval for a rate difference
based on Wald and Wilson confidence intervals for the individual rates are provided.}
\usage{confIntRateDiff(x, t, conf.level = 0.95)}
\arguments{
  \item{x}{vector of length 2, number of events in each group.}
  \item{t}{vector of length 2, total observation time in each group.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{rd}{Estimated rate difference.}
  \item{CIs}{Dataframe containing confidence intervals for the rate difference.}
}

\references{
Held, L., Rufibach, K. and Seifert, B. (2013).
Medizinische Statistik - Konzepte, Methoden, Anwendungen. 
Section 8.2. 
}

\author{Leonhard Held}
\seealso{\code{\link{wilsonRate}}}

\examples{
x <- c(30, 50)
t <- c(100, 120)
confIntRateDiff(x, t)$CIs
}

\keyword{htest}
