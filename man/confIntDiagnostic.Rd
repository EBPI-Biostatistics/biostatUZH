\name{confIntDiagnostic}
\alias{confIntDiagnostic}

\title{Compute confidence intervals for operating characteristics of a diagnostic test}
\description{
Compute confidence intervals for sensitivity, specificity, positive and
negative likelihood ratio and diagnostic odds ratio of a diagnostic test. Optionally also positive
and negative predictive value. 
}
\usage{confIntDiagnostic(tp, fp, tn, fn, conf.level = 0.95, cohort=FALSE, pr=NA, digits=NA)}
\arguments{
  \item{tp}{Number of true positives.}
  \item{fp}{Number of false positives.}
  \item{tn}{Number of true negatives.}
  \item{fn}{Number of false negatives.}
  \item{conf.level}{Confidence level for confidence intervals.}
  \item{cohort}{Logical indicator whether data comes from a cohort study.}
  \item{pr}{Prevalence.}
  \item{digits}{Number of digits.}
}

\value{
A dataframe containing the estimated confidence intervals for
sensitivity, specificity, positive and
negative likelihood ratio. Optionally also positive
and negative predictive value. }

\references{
Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\examples{
## Calculate confidence intervals for data from the Million Women Study

confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97)
confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, cohort=TRUE)
confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, pr=0.045)
confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, digits=2)

}


