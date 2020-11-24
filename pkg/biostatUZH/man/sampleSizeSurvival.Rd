\name{sampleSizeSurvival}
\alias{sampleSizeSurvival}
\alias{PrEvent}
\title{Compute sample size for a survival endpoint}
\description{Calculates two of the following quantities if one of them is given: 
  Sample size, required number of events and power.
  The distribution of survival times is assumed to be either non-parametric
  given a Kaplan-Meier estimate or following an exponential or Weibull model 
  given rate parameters or rate and shape parameters, respectively.
  The same time unit should be used for the survival functions, 
  the accrual period and the follow-up period.}
\usage{
sampleSizeSurvival(HR, a.length, f.length, sig.level=0.05, power=NULL,
                   n = NULL, n.events=NULL, alloc.ratio=1, drop.rate=0,
                   non.inf.margin=NULL, type="sup", dist="exp",
                   lambda=NULL, shape=NULL, survfit.ref=NULL,
                   alternative="two.sided", method="exact")
}
\arguments{
  \item{HR}{Hazard ratio of experimental group vs. control group.}
  \item{a.length}{Length of accrual period.}
  \item{f.length}{Length of follow-up period.}
  \item{sig.level}{Significance level.}
  \item{power}{Desired power.}
  \item{n}{Required sample size.}
  \item{n.events}{Required total number of events.}
  \item{alloc.ratio}{Allocation ratio: Ratio of the number of patients in the 
    experimental group divided by the number of patients in the control group.}
  \item{drop.rate}{Drop-out rate.}
  \item{non.inf.margin}{Non-inferiority margin.}
  \item{type}{Study type. Either "sup" or "superiority" for superiority
    studies or "non-inferiority" or "noninf" for non-inferiority
    studies. Default is "sup".}
  \item{dist}{Distributional assumption for survival times. Either "exp"
    or "exponential" for exponential, "weib" or "weibull" for Weibull,
    or "nonp" or "non-parametric" for non-parametric. Defaults to "exp".}
  \item{lambda}{The rate parameter in the exponential distribution and the scale
parameter in the Weibull distribution.}
  \item{shape}{The shape parameter of the Weibull distribution}
  \item{survfit.ref}{The survfit object, i.e. the Kaplan-Meier estimate
    of survival in the control group.}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
alternative hypothesis is one or two-sided.}
  \item{method}{In c("exact","approx","approximate"), depending on whether the integral for the
probability of event is solved with the function integrate
or approximated. If a non-parametric approach is chosen only
the approximate approach is available.}
}
\value{A list with the entries
  \item{n}{Required sample size.}
  \item{HR}{Hazard ratio of experimental group vs. control group.}
  \item{power}{Desired power.}
  \item{sig.level}{Significance level.}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
    alternative hypothesis is one- or two-sided.}
  \item{distribution}{Distributional assumption for survival times. Either "exp"
    or "exponential" for exponential, "weib" or "weibull" for Weibull,
    or "nonp" or "non-parametric" for non-parametric. }
  \item{PrEvent}{Probability that a patient will experience the event during the study.}
  \item{Events}{Required total number of events.}

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

## survival of females in lung cancer data
data(lung, package = "survival")
lung.female <- subset(lung, sex==2)
survObj <- Surv(time = lung.female$time, event=lung.female$status==2, type='right')
fit <- survfit(survObj ~ 1)
## exponential model
exp.reg <- survreg(survObj~1, dist = "exponential")
## Weibull model
weib.reg <- survreg(survObj~1, dist = "weibull")

## estimate the sample size with 5 different approaches

## Non-parametric
sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
                   f.length = 400, dist = "nonp", survfit.ref = fit,
                   type = "sup", alloc.ratio = 1, method = "approx")

## Exponential, approximate
exp.lambda <- 1/exp(coef(exp.reg))
sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
                   f.length = 400, dist = "exp",
                   lambda = exp.lambda, type = "sup",
                   alloc.ratio = 1, method = "approx")

## Exponential, exact
sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
                   f.length = 400, dist = "exp",
                   lambda = exp.lambda, type = "sup",
                   alloc.ratio = 1, method = "exact")

weib.scale <- unname(exp(coef(weib.reg)))
weib.shape <- unname(1/weib.reg$scale)
## Weibull, approximate
sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
                   f.length = 400, dist = "weib", lambda = weib.scale,
                   type = "sup", alloc.ratio = 1, shape = weib.shape,
                   method = "approx")

## Weibull, exact
sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
                   f.length = 400, dist = "weib", lambda = weib.scale,
                   type = "sup", alloc.ratio = 1, shape = weib.shape,
                   method = "exact")
}
\seealso{The function \code{\link{sampleSizeSurvival}} depends on
  \code{\link{NumEvents}} and \code{\link{PrEvent}}.}

