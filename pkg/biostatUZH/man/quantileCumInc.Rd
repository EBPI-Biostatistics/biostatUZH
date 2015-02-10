\name{quantileCumInc}
\alias{quantileCumInc}
\title{Compute arbitrary quantile of a Cumulative Incidence estimate}
\description{In survival analysis, often some quantiles (mainly the median) of an estimated cumulative incidence 
function are of interest. This function computes quantiles of cumulative incidence functions received
via \code{cuminc} in the package \pkg{cmprsk}.}
\usage{quantileCumInc(time, event, group, quant = 0.5)}
\arguments{
  \item{time}{Event times, censored or observed.}
  \item{event}{Event indicator.}
  \item{group}{Quantiles can be computed for several survival curves, defined by \code{group}.}
  \item{quant}{Quantile to be computed. Real number in \eqn{[0, 1]}.}
}
\value{Dataframe with quantiles estimates for each event-group combination.}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}

\seealso{\code{cuminc} and \code{timepoint} in the package \pkg{cmprsk}.}
\examples{
if (require("cmprsk")) {
  ## illustrate on simulated data
  set.seed(1977)
  n <- 180
  time <- rgamma(n, 2, 1)
  status <- factor(sample(rep(0:2, n / 3)), levels = 0:2, labels = 
      c("censored", "Time to 2nd tumor", "Death"))
  x <- factor(round(runif(n, 1, 2)), levels = 1:2, labels = 
      c("Tmt A", "Tmt B"))
  
  ## plot
  plot(cuminc(time, status, x))

  ## median time to event
  quantileCumInc(time, status, group = x, quant = 0.25)
}
}

\keyword{htest}
\keyword{survival}
