\name{sampleSizeMcNemar}
\alias{sampleSizeMcNemar}
\title{Compute sample size for McNemar test}
\description{Compute sample size to test the null hypothesis of the McNemar test.}
\usage{sampleSizeMcNemar(p1, p2, alpha = 0.05, power = 0.8)}
\arguments{
  \item{p1}{Assumed marginal probability of the first row.}
  \item{p2}{Assumed marginal probability of the first column.}
  \item{alpha}{Significance level.}
  \item{power}{Desired power, i.e. 1 - probability of type II error.}
}
\details{Given two samples of paired binary observations, e.g. results of some positive/negative 
test on the same experimental units we want to assess the null hypothesis whether the probability of positives by the first
method is equal that of positives by the second method.}
\value{
  \item{N}{Vector of min, max, and mid sample size as explained in Lachenbruch (1982).}
}
\references{
Lachenbruch (1992).
Sample Size for Studies based on McNemar's test
\emph{Stat. Med.}, \bold{11}, 1521--1525.
}
\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}

\examples{
# example from Lachenbruch (1982), Table II, first row
sampleSizeMcNemar(p1 = 0.8, p2 = 0.9, alpha = 0.05, power = 0.9)
}

\keyword{htest}

