\name{confIntIndependentAUCDiff}
\alias{confIntIndependentAUCDiff}

\title{Computes confidence interval for the difference in AUC based on two
  independent samples.}
\description{Computes confidence interval for the difference in the area
  under the curve based on two independent samples.}
\usage{confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type="Wald", conf.level=0.95)}
\arguments{
  \item{casesA}{Values of the continuous variable from Test A for the cases.}
  \item{controlsA}{Values of the continuous variable from Test A for the controls.}
  \item{casesB}{Values of the continuous variable from Test B for the cases.}
  \item{controlsB}{Values of the continuous variable from Test B for the controls.}
  \item{type}{"Wald" or "Logit". Defaults to "Wald".}
  \item{conf.level}{Confidence level for confidence interval.}
}
\details{For type="Wald", standard Wald confidence intervals are
  calculated for AUC of both tests and their difference. For
  type="logit", the substitution method is used based on the logit
  transformation for the AUC of both tests. The confidence interval for
  the difference in AUC is then calculated using Newcombe's method. }

\value{
  A matrix with estimate and confidence limits for AUC of the two tests
  and their difference.  
}
\references{
Newcombe, R.G. (1998).
Interval estimation for the difference between independent proportions: Comparison of eleven methods.
\emph{Stat. Med.}, \bold{17}, 873--890.

Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\examples{

set.seed(12345)
casesA <- rnorm(200, mean=2.5)
controlsA <- rnorm(100)
casesB <- rnorm(100, mean=1.5)
controlsB <- rnorm(200)

confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type="Wald")
confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type="Logit")

}

\keyword{univar}
\keyword{htest}
