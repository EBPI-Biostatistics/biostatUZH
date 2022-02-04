\name{confIntRate}
\alias{confIntRate}
\alias{waldRate}
\alias{wilsonRate}

\title{Compute confidence interval for a Poisson rate via different methods}
\description{
  Compute a confidence interval for a Poisson rate using several
  methods. The individual methods are also
  available as separate functions \code{waldRate} and \code{wilsonRate}.
}
\usage{confIntRate(x, t, conf.level = 0.95)}
\arguments{
  \item{x}{Number of events.}
  \item{t}{Total observation time.}
  \item{conf.level}{Confidence level for confidence interval.}
}

\value{
A list with the entries:
  \item{p}{Estimated rate.}
  \item{CIs}{Dataframe containing the estimated confidence intervals.}
}
\references{
Held, L., Rufibach, K. and Seifert, B. (2013).
Medizinische Statistik - Konzepte, Methoden, Anwendungen. 
Section 8.2. 
}

\author{Leonhard Held}
\seealso{Functions for some of the intervals provided here are available in \pkg{Hmisc} (see the examples).}

\examples{
## Calculate confidence bounds for a Poisson rate by different methods. 
x <- 1
t <- 3
ci <- confIntRate(x, t)$CIs
ci
p <- confIntRate(x, t)$rate

plot(0, 0, type = 'n', ylim = c(0, 3), xlim = c(0, 3), xlab = 'p', 
    ylab = '', yaxt = 'n')
lines(ci[1, 2:3], c(1, 1))
lines(ci[2, 2:3], c(2, 2))
points(p, 1, pch=19)
points(p, 2, pch=19)
text(0.5, 0.85, 'wald')
text(0.5, 1.85, 'wilson')
}
