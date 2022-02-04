\name{confIntProportion}
\alias{confIntProportion}
\alias{wald}
\alias{wilson}
\alias{agresti}
\alias{jeffrey}
\alias{jeffreys}
\alias{clopperPearson}

\title{Compute confidence interval for a binomial proportion via different methods}
\description{
  Compute a confidence interval for a binomial proportion using several
  asymptotic and exact methods. The individual methods are also
  available as separate functions \code{wald}, \code{wilson},
  \code{agresti}, \code{jeffreys}, and \code{clopperPearson}.
}
\usage{confIntProportion(x, n, conf.level = 0.95)}
\arguments{
  \item{x}{Number of successes.}
  \item{n}{Total number of trials.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{p}{Estimated proportion.}
  \item{CIs}{Dataframe containing the estimated confidence intervals.}
}
\references{
All the intervals provided in these functions are compared in:

Brown, L.D., Cai, T.T., DasGupta, A. (2001).
Interval Estimation for a Binomial Proportion.
\emph{Statistical Science}, \bold{16(2)}, 101--133.
}

\author{Kaspar Rufibach and Leonhard Held}
\seealso{Functions for some of the intervals provided here are available in \pkg{Hmisc} (see the examples).}

\examples{
## Calculate confidence bounds for a binomial parameter by different methods. 
x <- 50
n <- 100
ci <- confIntProportion(x, n)$CIs
ci

plot(0, 0, type = 'n', ylim = c(0, 7), xlim = c(0, 1), xlab = 'p', 
    ylab = '', yaxt = 'n')
lines(ci[1, 2:3], c(1, 1))
lines(ci[2, 2:3], c(2, 2))
lines(ci[3, 2:3], c(3, 3))
lines(ci[4, 2:3], c(4, 4))
lines(ci[5, 2:3], c(5, 5))
text(0.5, 0.85, 'wald')
text(0.5, 1.85, 'wilson')
text(0.5, 2.85, 'agresti')
text(0.5, 3.85, 'jeffreys')
text(0.5, 4.85, 'clopper')

## compare intervals to those received by the function binconf in Hmisc:
if (require("Hmisc")) {
    binconf(x, n, method = "asymptotic")   # Wald
    binconf(x, n, method = "wilson")       # Wilson
    binconf(x, n, method = "exact")        # Clopper-Pearson
}
}

\keyword{univar}
\keyword{htest}
