\name{sampleSizeSurvival}
\alias{sampleSizeSurvival}
\alias{PrEvent}
\title{Compute sample size for a survival endpoint}
\description{Function sampleSizeSurvival calculates up to two quantities
  out of the following: 
  sample size, required number of events and power.}
\usage{
sampleSizeSurvival(HR, a.length, f.length, sig.level=0.05, power=NULL,
                   n = NULL, n.events=NULL, alloc.ratio=1, drop.rate=0,
                   non.inf.margin=NULL, type="sup", dist="exp",
                   lambda=NULL, shape=NULL, survfit.ref=NULL,
                   alternative="two.sided", method="exact")
}
\arguments{
  \item{HR}{Hazard ratio.}
  \item{a.length}{Length of accrual period.}
  \item{f.length}{Length of follow-up period.}
  \item{sig.level}{Significance level.}
  \item{power}{Desired power, i.e. 1 - probability of type II error.}
  \item{n}{Required sample size.}
  \item{n.events}{Required number of events.}
  \item{alloc.ratio}{Allocation ratio: Number of allocated patients to
    the reference group relative to the treatment group.}
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
    of survival in the reference group.}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
alternative hypothesis is one or two-sided.}
  \item{method}{In c("exact","approx","approximate"), depending on whether the integral for the
probability of event is solved with the function integrate
or approximated. If a non-parametric approach is chosen only
the approximate approach is available.}
}
\value{A list with the entries
  \item{n}{Required sample size.}
  \item{HR}{Hazard ratio.}
  \item{power}{Desired power, i.e. 1 - probability of type II error.}
  \item{sig.level}{Significance level.}
  \item{alternative}{In c("one.sided","two.sided"), depending on whether the
alternative hypothesis is one or two-sided.}
  \item{distribution}{Distributional assumption for survival times. Either "exp"
    or "exponential" for exponential, "weib" or "weibull" for Weibull,
    or "nonp" or "non-parametric" for non-parametric. }
  \item{PrEvent}{The probability that a patient will not experience the event during the study.}
  \item{Events}{The required total number of events.}

}
\references{
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

