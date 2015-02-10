\name{formatPercent}
\alias{formatPercent}
\title{Format a Numeric Proportion}
\description{
  Takes a number and formats it as a percentage.
  Used by \code{\link{confIntIndependentProportion}}.
}

\usage{
formatPercent(x, digits = 1)
}

\arguments{
  \item{x}{numeric vector of proportions}
  \item{digits}{number of digits}
}

\author{Leonhard Held}

\examples{
formatPercent(c(0.115, 0.5))  # "11.5%" "50.0%"
}
