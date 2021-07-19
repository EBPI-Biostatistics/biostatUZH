## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function sampleSizeSurvival


#' Compute sample size for a survival endpoint
#' 
#' Calculates two of the following quantities if one of them is given: Sample
#' size, required number of events and power.  The distribution of survival
#' times is assumed to be either non-parametric given a Kaplan-Meier estimate
#' or following an exponential or Weibull model given rate parameters or scale
#' and shape parameters, respectively.  The same time unit should be used for
#' the survival times, the accrual period and the follow-up period.
#' 
#' 
#' @aliases sampleSizeSurvival PrEvent
#' @param HR Hazard ratio of experimental group vs. control group.
#' @param a.length Length of accrual period.
#' @param f.length Length of follow-up period.
#' @param sig.level Significance level.
#' @param power Desired power.
#' @param n Required sample size.
#' @param n.events Required total number of events.
#' @param alloc.ratio Allocation ratio: Ratio of the number of patients in the
#' experimental group divided by the number of patients in the control group.
#' @param drop.rate Drop-out rate.
#' @param non.inf.margin Non-inferiority margin.
#' @param type Study type. Either "sup" or "superiority" for superiority
#' studies or "noninf" or "non-inferiority" for non-inferiority studies.
#' Default is "sup".
#' @param dist Distributional assumption for survival times. Either "exp" or
#' "exponential" for exponential, "weib" or "weibull" for Weibull, or "nonp" or
#' "non-parametric" for non-parametric. Default is "exp".
#' @param lambda The rate parameter in the exponential distribution and the
#' scale parameter in the Weibull distribution.
#' @param shape The shape parameter of the Weibull distribution
#' @param survfit.ref The survfit object, i.e. the Kaplan-Meier estimate of
#' survival in the control group.
#' @param alternative In c("one.sided","two.sided"), depending on whether the
#' alternative hypothesis is one- or two-sided.
#' @param method In c("exact","approx","approximate"), depending on whether the
#' integral for the probability of event is solved with the function integrate
#' or approximated. If a non-parametric approach is chosen only the approximate
#' approach is available.
#' @return A list with the entries \item{n}{Required sample size.}
#' \item{HR}{Hazard ratio of experimental group vs. control group.}
#' \item{power}{Desired power.} \item{sig.level}{Significance level.}
#' \item{alternative}{In c("one.sided","two.sided"), depending on whether the
#' alternative hypothesis is one- or two-sided.}
#' \item{distribution}{Distributional assumption for survival times. Either
#' "exp" or "exponential" for exponential, "weib" or "weibull" for Weibull, or
#' "nonp" or "non-parametric" for non-parametric. } \item{PrEvent}{Probability
#' that a patient will experience the event during the study.}
#' \item{Events}{Required total number of events.}
#' @author Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@@uzh.ch}
#' @seealso The function \code{\link{sampleSizeSurvival}} depends on
#' \code{\link{NumEvents}} and \code{\link{PrEvent}}.
#' @references Collett, D. (2015). \emph{Modelling Survival Data in Medical
#' Research}.  New York: Chapman and Hall/CRC.
#' 
#' Schoenfeld, D.A. (1983). Sample-Size Formula for the Proportional-Hazards
#' Regression Model.  \emph{Biometrics}, \bold{39}, 499--503.
#' @examples
#' 
#' 
#' ## survival of females in lung cancer data
#' data(lung, package = "survival")
#' lung.female <- subset(lung, sex==2)
#' survObj <- Surv(time = lung.female$time, event=lung.female$status==2, type='right')
#' fit <- survfit(survObj ~ 1)
#' ## exponential model
#' exp.reg <- survreg(survObj~1, dist = "exponential")
#' ## Weibull model
#' weib.reg <- survreg(survObj~1, dist = "weibull")
#' 
#' ## estimate the sample size with 5 different approaches
#' 
#' ## Non-parametric
#' sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
#'                    f.length = 400, dist = "nonp", survfit.ref = fit,
#'                    type = "sup", alloc.ratio = 1, method = "approx")
#' 
#' ## Exponential, approximate
#' exp.lambda <- 1/exp(coef(exp.reg))
#' sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
#'                    f.length = 400, dist = "exp",
#'                    lambda = exp.lambda, type = "sup",
#'                    alloc.ratio = 1, method = "approx")
#' 
#' ## Exponential, exact
#' sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
#'                    f.length = 400, dist = "exp",
#'                    lambda = exp.lambda, type = "sup",
#'                    alloc.ratio = 1, method = "exact")
#' 
#' weib.scale <- unname(exp(coef(weib.reg)))
#' weib.shape <- unname(1/weib.reg$scale)
#' ## Weibull, approximate
#' sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
#'                    f.length = 400, dist = "weib", lambda = weib.scale,
#'                    type = "sup", alloc.ratio = 1, shape = weib.shape,
#'                    method = "approx")
#' 
#' ## Weibull, exact
#' sampleSizeSurvival(HR = .65, power=.90, a.length = 400,
#'                    f.length = 400, dist = "weib", lambda = weib.scale,
#'                    type = "sup", alloc.ratio = 1, shape = weib.shape,
#'                    method = "exact")
#' 
sampleSizeSurvival <- function(HR, a.length, f.length, sig.level=0.05, power=NULL,
                               n = NULL, n.events=NULL, alloc.ratio=1, drop.rate=0,
                               non.inf.margin=NULL, type="sup", dist="exp",
                               lambda=NULL, shape=NULL, survfit.ref=NULL,
                               alternative="two.sided", method="exact") {
    ## number of events or power needed (depending on given arguments)
    if(is.null(power) & is.null(n.events) & is.null(n)){
        stop(paste("either the power or the number of events",
                   "or the sample size need to be specified", sep = "\n"))
    }
    ## not allowed to overspecify number of events and power
    if((!is.null(power) & !is.null(n.events)) |
       (!is.null(power) & !is.null(n)) |
       (!is.null(n) & !is.null(n.events))){
        stop(paste("either the power or the number of events",
                   "or the sample size may be specified", sep = "\n"))
    }
    ## probability of event
    pr.event <- PrEvent(HR = HR, a.length = a.length, f.length = f.length,
                        dist = dist, lambda = lambda, shape = shape, method = method,
                        survfit.ref = survfit.ref, alloc.ratio = alloc.ratio)
    ## power calculation
    if(is.null(power)){
        if(is.null(n.events)){
            n.events <- ceiling(pr.event * n)
        }
        power <- NumEvents(HR = HR, sig.level = sig.level, power = power,
                          n.events = n.events, alloc.ratio = alloc.ratio,
                          non.inf.margin = non.inf.margin, type = type,
                          alternative = alternative)
    }
    ## sample size calculation
    if(is.null(n)){
        if(is.null(n.events)){
            n.events <- NumEvents(HR = HR, sig.level = sig.level, power = power,
                                 n.events = n.events, alloc.ratio = alloc.ratio,
                                 non.inf.margin = non.inf.margin, type = type,
                                 alternative = alternative)
        }
        n <- n.events/pr.event
        n <- ceiling(n/(1-drop.rate))
    }
    ## return statement
    str <- structure(list(n=n, HR = HR, power=power, sig.level=sig.level,
                          alternative=alternative, distribution = dist,
                          PrEvent = pr.event, Events=n.events,
                          NOTE="n is the total sample size"))
    return(str)
}
