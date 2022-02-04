\name{confIntIndependentDiagnostic}
\alias{confIntIndependentDiagnostic}

\title{Compute confidence intervals for the comparison of two diagnostic
  tests from unpaired data}
\description{
  Compute confidence intervals for relative fractions and relative
  likelihood ratios. Specifically, confidence intervals for the relative true positive, true
  negative, false positive and false negative fraction as well as the
  relative positive and negative likelihood ratio are provided. 
}
\usage{confIntIndependentDiagnostic(tp, fp, tn, fn, conf.level = 0.95, adjust=FALSE)}
\arguments{
  \item{tp}{Number of true positives of the two tests.}
  \item{fp}{Number of false positives.}
  \item{tn}{Number of true negatives.}
  \item{fn}{Number of false negatives.}
  \item{conf.level}{Confidence level for confidence intervals.}
  \item{adjust}{Adjusted CIs? by default FALSE.}
}

\value{
A dataframe containing the estimated confidence intervals for
the measures of relative accuracy (first versus second test).}

\references{
Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\examples{
## Calculate confidence intervals for data from a (hypothetical)
## randomized unpaired study of chorionic villus sampling (CVS)
## versus early amniocentesis (EA) for fetal abnormality from
## Pepe (2003)

tp <- c(116, 111)
fp <- c(34, 111)
tn <- c(4844, 4765)
fn <- c(6, 13)

print(confIntIndependentDiagnostic(tp=tp, fp=fp, tn=tn, fn=fn))

}


