\name{standardErrorAUCDiff}
\alias{standardErrorAUCDiff}

\title{Standard Error of a AUC difference}
\description{
  Computes the standard error of the difference in area under the curve for paired samples. 
  
}
\usage{standardErrorAUCDiff(cases, controls)}
\arguments{
  \item{cases}{Matrix with values of the continuous variable for the
    cases. First column gives values of test A, second gives values of
    test B.}
  \item{controls}{Matrix with values of the continuous variable for the controls. First column gives values of test A, second gives values of
    test B.}
}

\value{
The standard error. 
}
\references{
Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Leonhard Held}

\keyword{univar}
\keyword{htest}
