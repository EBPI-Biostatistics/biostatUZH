\name{formatPval}
\alias{formatPval}

\title{Format p-Values}

\description{Formats p-values in a desired manner.}

\usage{formatPval(x, break.eps = 1e-04, break.middle = 0.01, ...)}

\arguments{
   \item{x}{a numeric vector.}
  \item{break.eps}{numeric value, see section details.}
  \item{break.middle}{numeric value, see section details.}
  \item{...}{additional arguments passed to \code{\link{format}}.}
}

\value{
A vector of \code{length(x)} with formatted p-values. 
}

\details{
  \itemize{
  \item p-values \code{< break.eps} are written as: \code{"< break.eps"}.
  \item p-values \code{>= break.eps} but \code{< break.middle} have 1 digit.
  \item p-values \code{>= break.middle} have 2 digits.
  }
}

\author{Sina Rueeger and Sebastian Meyer}

\seealso{
  the \pkg{base} function \code{\link{format.pval}},
  \code{\link[Hmisc]{format.pval}} in package \pkg{Hmisc},
  \code{\link[reporttools]{formatPval}} in package \pkg{reporttools},
  \code{\link[surveillance]{formatPval}} in package \pkg{surveillance}
}

\examples{
x <- c(1e-8, 0.00568, 0.0345, 0.885)
biostatUZH::formatPval(x)  # "< 0.0001" "0.006" "0.035" "0.89"

## compare to formatting of other packages
reporttools::formatPval(x) # "< 0.0001" "0.0057" "0.03" "0.88"
if (requireNamespace("Hmisc")) {
    Hmisc::format.pval(x)  # "0" "0.00568" "0.03450" "0.88500"
}
if (requireNamespace("surveillance")) {
    surveillance::formatPval(x)  # "<0.0001" "0.0057" "0.035" "0.89"
}

## adapt break.middle
biostatUZH::formatPval(x, break.middle = 0.001)
}
