\name{confIntAUC}
\alias{confIntAUC}

\title{Compute confidence interval for AUC}
\description{
  Compute Wald confidence intervals for the area under the curve on the
  original and on the logit scale. 
  
}
\usage{confIntAUC(cases, controls, conf.level = 0.95)}
\arguments{
  \item{cases}{Values of the continuous variable for the cases.}
  \item{controls}{Values of the continuous variable for the controls.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{p}{Estimated proportion.}
  \item{CIs}{Dataframe containing the estimated confidence intervals.}
}
\references{

The computations follow Chapter 10 in

Altman, D.G., Machin, D., Bryant, T.N. and Gardner, M.J. (2001).
\emph{Statistics with confidence}. 2nd Edition, 2000. BMJ Books.  
}

\author{Leonhard Held}
\seealso{Earlier function with this name is now called summaryROC().}

\examples{
set.seed(12345)
cases <- rnorm(100, mean=2)
controls <- rnorm(50)
(confIntAUC(cases, controls))
}

\keyword{univar}
\keyword{htest}
