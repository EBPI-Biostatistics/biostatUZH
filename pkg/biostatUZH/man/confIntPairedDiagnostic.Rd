\name{confIntPairedDiagnostic}
\alias{confIntPairedDiagnostic}

\title{Compute confidence intervals for the comparison of two diagnostic
  tests from paired data}
\description{
  Compute confidence intervals for relative fractions. Specifically,
  confidence intervals 
  for the relative true positive, true
  negative, false positive and false negative fraction are provided. 
}
\usage{confIntPairedDiagnostic(Diseased, nonDiseased, conf.level = 0.95, adjust=FALSE)}
\arguments{
  \item{Diseased}{Frequency table with results from both tests in the
  diseased population. First row and first column refers to negative
  test results.}
  \item{nonDiseased}{Frequency table with results from both tests in the non-diseased population. First row and first column refers to negative
  test results.}
  \item{conf.level}{Confidence level for confidence intervals.}
  \item{adjust}{Adjusted CIs? By default FALSE.}
}

\value{
A dataframe containing the estimated confidence intervals for
the measures of relative accuracy (column test vs. row test).}

\references{
Section 3.3 in Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\examples{
## Calculate confidence intervals for data from a study performing the 
## Exercise stress test (EST) and the determination of chest pain
## history (CPH) for patients with suspected or probably coronary
## heart disease (CHD).

Diseased <- matrix(c(25, 183, 29, 786), ncol=2, nrow=2, byrow=TRUE)
nonDiseased <- matrix(c(151, 176, 46, 69), ncol=2, nrow=2, byrow=TRUE)
colnames(Diseased) <- colnames(nonDiseased) <- c("CPH=0", "CPH=1")
rownames(Diseased) <- rownames(nonDiseased) <- c("EST=0", "EST=1")

print(confIntPairedDiagnostic(Diseased=Diseased, nonDiseased=nonDiseased))

}


