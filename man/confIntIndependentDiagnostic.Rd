% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confIntDiagnostic.R
\name{confIntIndependentDiagnostic}
\alias{confIntIndependentDiagnostic}
\title{Confidence intervals for the comparison of two diagnostic tests from
unpaired data}
\usage{
confIntIndependentDiagnostic(
  tp,
  fp,
  tn,
  fn,
  names = c(1, 2),
  conf.level = 0.95,
  adjust = FALSE
)
}
\arguments{
\item{tp}{Vector of length 2. Number of true positives of the two tests.}

\item{fp}{Vector of length 2. Number of false negatives of the two tests.}

\item{tn}{Vector of length 2. Number of true negatives of the two tests.}

\item{fn}{Vector of length 2. Number of false negatives of the two tests.}

\item{names}{Vector of length 2. Names of first and second diagnostic tests.
Default to c(1, 2).}

\item{conf.level}{Confidence level for confidence interval. Default is 0.95.}

\item{adjust}{Logical of length one indicating whether to compute adjusted
CIs or not. Default is FALSE.}
}
\value{
A data.frame containing the estimated sensitivity, specificity,
LR+ and LR- in the two tests, as well as the ratio (first versus second test)
with the corresponding confidence interval.
}
\description{
Compute confidence intervals for relative sensitivity, relative specificity and
relative likelihood ratios.
}
\examples{

## Calculate confidence intervals for data from a (hypothetical)
## randomized unpaired study of early amniocentesis (EA)
## versus chorionic villus sampling (CVS)
## for fetal abnormality from
## Pepe (2003)

tp <- c(116, 111)
fp <- c(34, 111)
tn <- c(4844, 4765)
fn <- c(6, 13)
confIntIndependentDiagnostic(tp = tp, fp = fp, tn = tn, fn = fn,
names = c("EA", "CVS"))

}
\references{
Pepe, M.S. (2003) \emph{The statistical evaluation of medical
tests for classification and prediction}. Oxford: Oxford University Press.
}
\author{
Leonhard Held, Charlotte Micheloud
}
