\name{printResults}
\alias{printResults}
\title{Confidence interval for a given effect}
\description{Prints the confidence interval and the p-value for a given effect and standard error.}

\usage{
printResults(theta, se.theta, conf.level=0.95, digits=3)
}

\arguments{
  \item{theta}{the effect size.}
  \item{se.theta}{the standard error.}
  \item{conf.level}{Confidence level. Defualt is 95\%.}
  \item{digits}{number of digits for the effect.}
}

\value{
  Prints effect and associated confidence p-values.
}

\author{Leonhard Held}

\examples{
printResults(6.7, 2.1)
}
