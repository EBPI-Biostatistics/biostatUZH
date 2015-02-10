\name{confIntPairedProportion}
\alias{confIntPairedProportion}
\alias{McNemar}
\alias{paired}

\title{Compute confidence interval for the difference of paired binomial proportions using Newcombe's method}
\description{Compute confidence interval for the difference of paired binomial proportions using Newcombe's method.}
\usage{confIntPairedProportion(x, conf.level = 0.95)}
\arguments{
  \item{x}{A two-dimensional contingency table in matrix form.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{p1}{Estimated proportion \eqn{p_{1+}}.}
  \item{p2}{Estimated proportion \eqn{p_{+1}}.}
  \item{newcombeCI}{Confidence interval for the difference of paired proportions, computed according to Newcombe (1998).}
}
\references{
The Newcombe interval is introduced in 

Newcombe, R.G. (1998).
Improved confidence intervals for the difference between binomial proportions based on paired data.
\emph{Stat. Med.}, \bold{17}, 2635--2650.

A worked out example can be found in 

Altman, D.G., Machin, D., Bryant, T.N., Gardner, M.J. (2000).
Statistics with confidence.
University Press Belfast.
}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}

\examples{
# Calculate confidence interval for the example in Altman et al (2000), Table 6.2
altman62 <- rbind(c(14, 5), c(0, 22))
confIntPairedProportion(x = altman62)

\dontrun{
# exact McNemar test
library(exact2x2)
mcnemar.exact(altman62)
mcnemar.test(altman62)
}
}

\keyword{univar}
\keyword{htest}
