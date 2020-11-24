\name{NumEvents}
\alias{NumEvents}
\title{Compute number of events for a survival endpoint}
\description{Calculates either the number of events or the power.}
\usage{
NumEvents(HR, sig.level=0.05, power=NULL, n.events=NULL, alloc.ratio=1,
non.inf.margin=NULL, type="sup", alternative="two.sided")
}
\arguments{
  \item{HR}{Hazard ratio.}
  \item{sig.level}{Significance level.}
  \item{power}{Desired power.}
  \item{n.events}{Required number of events.}
  \item{alloc.ratio}{Allocation ratio: Ratio of the number of patients in the 
    experimental group divided by the number of patients in the control group.}
  \item{non.inf.margin}{Non-inferiority margin.}
  \item{type}{Study type. Either "sup" or "superiority" for superiority
    studies or "non-inferiority" or "noninf" for non-inferiority
    studies. Default is "sup".}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
    alternative hypothesis is one- or two-sided.}
}
\value{Returns either the required total number of events or the power. 
  }
\references{
Collett, D. (2015).
\emph{Modelling Survival Data in Medical Research}. 
New York: Chapman and Hall/CRC.

Schoenfeld, D.A. (1983).
Sample-Size Formula for the Proportional-Hazards Regression Model. 
\emph{Biometrics}, \bold{39}, 499--503.
}
\author{Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@uzh.ch}}

\examples{
NumEvents(HR = 0.65, sig.level = 0.05, power = 0.9,
alloc.ratio = 1, type = "sup", alternative = "two.sided")
}


