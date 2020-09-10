\name{NumEvents}
\alias{NumEvents}
\title{Compute number of events for a survival endpoint}
\description{Function NumEvents calculates either the number of events or
  the power.}
\usage{
NumEvents(HR, sig.level=0.05, power=NULL, n.events=NULL, alloc.ratio=1,
non.inf.margin=NULL, type="sup", alternative="two.sided")
}
\arguments{
  \item{HR}{Hazard ratio.}
  \item{sig.level}{Significance level.}
  \item{power}{Desired power, i.e. 1 - probability of type II error.}
  \item{n.events}{Required number of events.}
  \item{alloc.ratio}{Allocation ratio: Number of allocated patients to
    the reference group relative to the treatment group.}
  \item{non.inf.margin}{Non-inferiority margin.}
  \item{type}{Study type. Either "sup" or "superiority" for superiority
    studies or "non-inferiority" or "noninf" for non-inferiority
    studies. Default is "sup".}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
alternative hypothesis is one or two-sided.}
}
\value{Returns either the required total number of events or the power. 
  }
\references{
}
\author{Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@uzh.ch}}

\examples{
NumEvents(HR = 0.65, sig.level = 0.05, power = 0.9,
alloc.ratio = 1, type = "sup", alternative = "two.sided")
}


