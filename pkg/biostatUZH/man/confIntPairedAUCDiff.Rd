\name{confIntPairedAUCDiff}
\alias{confIntPairedAUCDiff}

\title{Computes confidence interval for the difference in AUC based on
  two paired samples.}

\description{ Computes confidence interval for the
difference in the area under the curve based on two paired samples.
}

\usage{confIntAUC(cases, controls, conf.level = 0.95)}
\arguments{
  \item{cases}{Matrix with values of the continuous variable for the
    cases. First column gives values of test A, second gives values of
    test B.}
  \item{controls}{Matrix with values of the continuous variable for the controls. First column gives values of test A, second gives values of
    test B.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
  A matrix with estimate and confidence limits for AUC of the two tests
  and their difference.  
}
\references{
Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\examples{

data(wiedat2b)
ind <- wiedat2b[,3]
cases <- wiedat2b[(ind==1),1:2]
controls <- wiedat2b[(ind==0),1:2]

confIntPairedAUCDiff(cases, controls)

}

\keyword{univar}
\keyword{htest}
