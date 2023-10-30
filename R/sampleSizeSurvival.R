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
#' @param dist Distributional assumption for survival times. Either "exponential",
#' "weibull" or "non-parametric". Default is "exponential".
#' @param lambda The rate parameter in the exponential distribution and the
#' scale parameter in the Weibull distribution.
#' @param shape The shape parameter of the Weibull distribution
#' @param survfit.ref The survfit object, i.e. the Kaplan-Meier estimate of
#' survival in the control group.
#' @param alternative Either "one.sided" or "two.sided" depending on whether the
#' alternative hypothesis is one- or two-sided.
#' @param method Either "exact" or "approx", indicating whether the
#' integral for the probability of event is solved by integration the function
#' or approximated. If a non-parametric approach is chosen only the approximate
#' approach is available.
#' @return A list with the entries
#' \item{n}{Required total sample size.}
#' \item{HR}{Hazard ratio of experimental group vs. control group.}
#' \item{power}{Power of the study.}
#' \item{sig.level}{Significance level.}
#' \item{alternative}{Indicates whether the the alternative hypothesis is one- or two-sided.}
#' \item{distribution}{Distributional assumption for survival times.
#' Either "exponential" for exponential, "weibull" or "non-parametric".}
#' \item{PrEvent}{Probability that a patient will experience the event during the study.}
#' \item{Events}{Required total number of events.}
#' @author Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@@uzh.ch}
#' @seealso The function \code{\link{sampleSizeSurvival}} depends on
#' \code{\link{survEvents}} and \code{\link{PrEvent}}.
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
#' lung.female <- subset(lung, sex == 2)
#' survObj <- Surv(time = lung.female$time, event=lung.female$status == 2,
#'                 type='right')
#' fit <- survfit(survObj ~ 1)
#' ## exponential model
#' exp.reg <- survreg(survObj~1, dist = "exponential")
#' ## Weibull model
#' weib.reg <- survreg(survObj~1, dist = "weibull")
#'
#' ## estimate the sample size with 5 different approaches
#'
#' ## Non-parametric
#' sampleSizeSurvival(HR = 0.65, power = 0.9, a.length = 400,
#'                    f.length = 400, dist = "non-parametric", survfit.ref = fit,
#'                    alloc.ratio = 1, method = "approx")
#'
#' ## Exponential, approximate
#' exp.lambda <- 1/exp(coef(exp.reg))
#' sampleSizeSurvival(HR = 0.65, power = 0.90, a.length = 400,
#'                    f.length = 400, dist = "exponential",
#'                    lambda = exp.lambda,
#'                    alloc.ratio = 1, method = "approx")
#'
#' ## Exponential, exact
#' sampleSizeSurvival(HR = 0.65, power = 0.90, a.length = 400,
#'                    f.length = 400, dist = "exponential",
#'                    lambda = exp.lambda,
#'                    alloc.ratio = 1, method = "exact")
#'
#' weib.scale <- unname(exp(coef(weib.reg)))
#' weib.shape <- unname(1/weib.reg$scale)
#' ## Weibull, approximate
#' sampleSizeSurvival(HR = 0.65, power = 0.90, a.length = 400,
#'                    f.length = 400, dist = "weibull", lambda = weib.scale,
#'                    alloc.ratio = 1, shape = weib.shape,
#'                    method = "approx")
#'
#' ## Weibull, exact
#' sampleSizeSurvival(HR = 0.65, power = 0.90, a.length = 400,
#'                    f.length = 400, dist = "weibull", lambda = weib.scale,
#'                    alloc.ratio = 1, shape = weib.shape,
#'                    method = "exact")
#'
#' @export
sampleSizeSurvival <- function(
    HR, a.length, f.length, sig.level = 0.05, power = NULL,
    n = NULL, n.events = NULL, alloc.ratio = 1, drop.rate = 0,
    non.inf.margin = 0,
    dist = c("exponential", "weibull", "non-parametric"),
    lambda = NULL, shape = NULL, survfit.ref = NULL,
    alternative = c("two.sided", "one.sided"), method = c("exact", "approx")
) {

    stopifnot(!is.null(dist))
    dist <- match.arg(dist)
    stopifnot(!is.null(alternative))
    alternative <- match.arg(alternative)
    stopifnot(!is.null(method))
    method <- match.arg(method)

    ## number of events or power needed (depending on given arguments)
    if (is.null(power) && is.null(n.events) && is.null(n)) {
        stop(paste("either the power or the number of events",
                   "or the sample size need to be specified", sep = "\n"))
    }
    ## not allowed to overspecify number of events and power
    if ((!is.null(power) && !is.null(n.events)) ||
       (!is.null(power) && !is.null(n)) ||
       (!is.null(n) && !is.null(n.events))) {
        stop(paste("either the power or the number of events",
                   "or the sample size may be specified", sep = "\n"))
    }

    ## more checks on arguments are done in PrEvent and survEvents


    ## probability of event
    pr.event <- PrEvent(HR = HR, a.length = a.length, f.length = f.length,
                        dist = dist, lambda = lambda, shape = shape, method = method,
                        survfit.ref = survfit.ref, alloc.ratio = alloc.ratio)
    ## power calculation
    if (is.null(power)) {
        if (is.null(n.events)) {
            n.events <- ceiling(pr.event * n)
        }
        power <- survEvents(HR = HR, sig.level = sig.level, power = power,
                            n.events = n.events, alloc.ratio = alloc.ratio,
                            non.inf.margin = non.inf.margin,
                            alternative = alternative)
    }
    ## sample size calculation
    if (is.null(n)) {
        if (is.null(n.events)) {
            n.events <- survEvents(
                HR = HR, sig.level = sig.level, power = power,
                n.events = n.events,
                alloc.ratio = alloc.ratio,
                non.inf.margin = non.inf.margin,
                alternative = alternative
            )
        }
        n <- n.events / pr.event
        n <- ceiling(n / (1 - drop.rate))
    }

    list(n = n, HR = HR, power = power, sig.level = sig.level,
         alternative = alternative, distribution = dist,
         PrEvent = pr.event, Events = n.events)
}

#' @importFrom utils tail
PrEvent <- function(
    HR, a.length, f.length,
    dist = c("exponential", "weibull", "non-parametric"),
    lambda = NULL, shape = NULL, survfit.ref = NULL, alloc.ratio = 1,
    method = c("exact", "approx")
) {

    stopifnot(is.numeric(HR), length(HR) > 0,
              is.finite(HR),
              0 < HR,

              !is.null(dist))
    dist <- match.arg(dist)
    stopifnot(!is.null(method))
    method <- match.arg(method)

    p1 <- alloc.ratio / (alloc.ratio + 1)
    p2 <- 1 - p1
    a <- a.length
    f <- f.length
    if (dist == "exponential") {
        stopifnot(!is.null(lambda), lambda > 0)
        if (method == "exact") {
            pr.event.cnt <- 1 - 1 / a  *
                integrate(function(x) {
                    1 - pexp(x, lambda)
                }, lower = f, upper = a + f)$value
            pr.event.trt <- 1 - 1 / a  *  integrate(function(x) {
                (1 - pexp(x, lambda))^HR
            }, lower = f, upper = a + f)$value
        } else {
            S1 <- 1 - pexp(f, lambda)
            S2 <- 1 - pexp(.5 * a + f, lambda)
            S3 <- 1 - pexp(a + f, lambda)
        }
    }
    if (dist == "weibull") {
        stopifnot(!is.null(lambda), !is.null(shape))
        if (method == "exact") {
            pr.event.cnt <- 1 - 1 / a * integrate(function(x) {
                1 - pweibull(x, shape, lambda)
            }, lower = f, upper = a + f)$value
            pr.event.trt <- 1 - 1 / a * integrate(function(x) {
                (1 - pweibull(x, shape, lambda))^HR
            }, lower = f, upper = a + f)$value
        } else {
            S1 <- 1 - pweibull(f, shape, lambda)
            S2 <- 1 - pweibull(.5 * a + f, shape, lambda)
            S3 <- 1 - pweibull(a + f, shape, lambda)
        }
    }
    if (dist == "non-parametric") {
        stopifnot(!is.null(survfit.ref))
        if (method == "exact") {
            method <- "approx"
        }
        time <- survfit.ref$time
        surv <- survfit.ref$surv
        index1 <- utils::tail(which(f >= time), 1)
        index2 <- utils::tail(which(.5 * a + f >= time), 1)
        index3 <- utils::tail(which(a + f >= time), 1)
        if (is.na(index3)) {
            stop("accrual + follow-up length longer than duration of survfit.ref")
        }
        S1 <- surv[index1]
        S2 <- surv[index2]
        S3 <- surv[index3]
    }
    if (method == "approx") {
        pr.event.cnt <- 1 - 1 / 6 * (S1 + 4 * S2 + S3)
        pr.event.trt <- 1 - 1 / 6 * (S1^HR + 4 * S2^HR + S3^HR)
    }
    pr.event <- p1 * pr.event.trt + p2 * pr.event.cnt
    return(pr.event)
}
