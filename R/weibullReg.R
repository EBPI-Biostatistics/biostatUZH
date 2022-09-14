#' Weibull regression for survival data
#' 
#' \code{weibullReg} performs Weibull regression using the
#' \code{\link{survreg}} function, and transforms the estimates to a more
#' natural parameterization. Additionally, it produces hazard ratios
#' (corresponding to the proportional hazards interpretation), and event time
#' ratios (corresponding to the accelerated failure time interpretation) for
#' all covariates.
#' 
#' Details regarding the transformations of the parameters and their standard
#' errors can be found in Klein and Moeschberger (2003, Chapter 12). An
#' explanation of event time ratios for the accelerated failure time
#' interpretation of the model can be found in Carroll (2003). A general
#' overview can be found in the \code{vignette("weibull")} of this package, or
#' in the documentation for \code{\link{survreg2weibull}}.
#' 
#' @param formula A \code{\link{Surv}} formula.
#' @param data The dataset containing all variables referenced in
#' \code{formula}.
#' @param conf.level Specifies that \eqn{1-\alpha} level confidence intervals
#' for the hazard and event time ratios should be produced.
#' @return \item{formula}{The formula for the Weibull regression model.}
#' \item{coef}{The transformed maximum likelihood estimates, with standard
#' errors.} \item{HR}{The hazard ratios for each of the predictors, with
#' \eqn{1-\alpha} level confidence intervals.} \item{ETR}{The event time ratios
#' (acceleration factors) for each of the predictors, with \eqn{1-\alpha} level
#' confidence intervals.} \item{summary}{The summary output from the original
#' \code{\link{survreg}} model.}
#' @author Sarah R. Haile
#' @seealso \code{\link{survreg2weibull}}, \code{\link[survival]{survreg}}.
#' @references Carroll, K. (2003).  On the use and utility of the Weibull model
#' in the analysis of survival data. \emph{Controlled Clinical Trials},
#' \bold{24}, 682--701.
#' 
#' Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis: techniques
#' for censored and truncated data}.  Springer.
#' @examples
#' 
#' data(larynx)
#' weibullReg(Surv(time, death) ~ factor(stage) + age, data = larynx)
#' 
#' @import survival
#' @export
weibullReg <- function (formula, data, conf.level = 0.95) 
{
    ## 'formula' checked in survreg()
    ## 'data' and 'conf.level' checked in survreg2weibull()
    m <- survreg(formula = formula, data = data, dist = "weibull")
    mle <- survreg2weibull(model = m, conf.level = conf.level)
    list(formula = formula, coef = mle$vars, HR = mle$HR, 
         ETR = mle$ETR, summary = summary(m))
}


#' Transformation of survreg output for the Weibull distribution
#' 
#' Transforms output from \code{\link{survreg}} using the Weibull distribution
#' to a more natural parameterization. See details for more information.
#' 
#' The \code{\link{survreg}} function fits a Weibull accelerated failure time
#' model of the form
#' 
#' \deqn{\log t = \mu + \gamma^T Z + \sigma W,}
#' 
#' where \eqn{Z} is a matrix of covariates, and \eqn{W} has the extreme value
#' distribution, \eqn{\mu} is the intercept, \eqn{\gamma} is a vector of
#' parameters for each of the covariates, and \eqn{\sigma} is the scale. The
#' usual parameterization of the model, however, is defined by hazard function
#' 
#' \deqn{h(t|Z) = \alpha \lambda t^{\alpha - 1} \exp(\beta^T Z).}
#' 
#' The transformation is as follows: \eqn{\alpha = 1/\sigma}, \eqn{\lambda =
#' \exp(-\mu/\sigma)}, and \eqn{\beta=-\gamma/\sigma}, and estimates of the
#' standard errors can be found using the delta method.
#' 
#' The Weibull distribution has the advantage of having two separate
#' interpretations. The first, via proportional hazards, leads to a hazard
#' ratio, defined by \eqn{\exp \beta}. The second, of accelerated failure
#' times, leads to an event time ratio (also known as an acceleration factor),
#' defined by \eqn{\exp (-\beta/\alpha)}.
#' 
#' Further details regarding the transformations of the parameters and their
#' standard errors can be found in Klein and Moeschberger (2003, Chapter 12).
#' An explanation of event time ratios for the accelerated failure time
#' interpretation of the model can be found in Carroll (2003). A general
#' overview can be found in the \code{vignette("weibull")} of this package.
#' 
#' @param model A \code{\link{survreg}} model, with \code{dist = "weibull"}
#' (the default).
#' @param conf.level Significance level used to produce two-sided
#' \eqn{1-\alpha/2} confidence intervals for the hazard and event time ratios.
#' @return \item{vars}{A matrix containing the values of the transformed
#' parameters and their standard errors} \item{HR}{A matrix containing the
#' hazard ratios for the covariates, and \eqn{1-\code{level}/2} confidence
#' intervals.} \item{ETR}{A matrix containing the event time ratios for the
#' covariates, and \eqn{1-\code{conf.level}/2} confidence intervals.}
#' @author Sarah R. Haile
#' @seealso Requires the packages \pkg{survival} and \pkg{prodlim}. This
#' function is used by \code{\link{weibullReg}}.
#' @references Carroll, K. (2003).  On the use and utility of the Weibull model
#' in the analysis of survival data. \emph{Controlled Clinical Trials},
#' \bold{24}, 682--701.
#' 
#' Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis: techniques
#' for censored and truncated data}. Springer.
#' @keywords survival regression
#' @aliases ConvertWeibull
#' @examples
#' 
#' data(larynx)
#' survreg2weibull(survreg(Surv(time, death) ~ stage + age, data = larynx),
#'                 conf.level = 0.95)
#' 
#' @export
survreg2weibull <- function (model, conf.level = 0.95)
{
    stopifnot(inherits(x = model, what = "survreg"),
              !is.null(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1)
    alpha <- 1 - conf.level
    qa <- qnorm(1 - alpha / 2)
    Int.Only <- (nrow(summary(model)$table) == 2)
    sigma <- summary(model)$scale
    mu <- summary(model)$coef[1]
    k <- length(summary(model)$coef) - 1
    if (!Int.Only) {
        gamma <- summary(model)$coef[2:(k + 1)]
    }
    lambda <- exp(-mu / sigma)
    alpha <- 1 / sigma
    tmp <- c(lambda, alpha)
    names(tmp) <- c("lambda", "alpha")
    if (!Int.Only) {
        beta <- -gamma / sigma
        tmp <- c(lambda, alpha, beta)
        names(tmp) <- c("lambda", "alpha",
                        names(summary(model)$coef[2:(k + 1)]))
    }
    var1 <- summary(model)$var
    var.mu <- diag(var1)[1]
    var.sigma <- var1[(k + 2), (k + 2)] * exp(2 * log(sigma))
    if (!Int.Only) {
        var.gamma <- var1[2:(k + 1), 2:(k + 1)]
        if(k > 1)
            var.gamma <- diag(var.gamma)
        se.gamma <- sqrt(var.gamma)
    }
    cov.mu.sigma <- var1[(k + 2), 1] * sigma
    var.alpha <- var.sigma / (sigma^4)
    var.lambda <- exp(-2 * mu / sigma) *
        (var.mu / sigma^2 - 2 * mu / sigma^3 * cov.mu.sigma +
         mu^2 / sigma^4 * var.sigma)
    var <- c(sqrt(var.lambda), sqrt(var.alpha))
    if (!Int.Only) {
        cov.gamma.sigma <- var1[2:(k + 1), (k + 2)] * sigma
        var.beta <- 1 / sigma^2 * (var.gamma - (2 * gamma / sigma) * 
            (cov.gamma.sigma) + (((gamma / sigma)^2) * var.sigma))
        se.beta <- sqrt(var.beta);
        var <- c(sqrt(var.lambda), sqrt(var.alpha), se.beta)
        HR <- cbind(HR = exp(beta), LB = exp(beta - qa * se.beta), 
            UB = exp(beta + qa * se.beta))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
        ETR <- cbind(ETR = exp(gamma),
                     LB = exp(gamma - qa * se.gamma),
                     UB = exp(gamma + qa * se.gamma))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
    }
    tmp1 <- rbind(tmp, var)
    rownames(tmp1) <- c("Estimate", "SE")
    ret <- list(vars = t(tmp1))
    if (!Int.Only) {
        ret$HR = HR
        ret$ETR = ETR
    }
    ret
}



#' @export
ConvertWeibull <- function(model, conf.level){
    .Deprecated(new = "survreg2weibull")
    survreg2weibull(model = model, conf.level = conf.level)
}



#' Diagnostic plot of adequacy of Weibull distribution
#' 
#' This function constructs a diagnostic plot of the adequacy of the Weibull
#' distribution for survival data with respect to one categorical covariate. If
#' the Weibull distribution fits the data well, then the lines produced should
#' be linear and parallel.
#' 
#' As discussed in Klein and Moeschberger (2003), one method for checking the
#' adequacy of the Weibull model with a categorical covariate is to produce
#' stratified Kaplan-Meier estimates (KM), which can be transformed to estimate
#' the log cumulative hazard for each stratum. Then in a plot of \eqn{\log(t)}
#' versus \eqn{\log(-\log(KM))}, the lines should be linear and parallel. This
#' can be seen as the log cumulative hazard for the Weibull distribution is
#' 
#' \deqn{\log H(t) = \log \lambda + \alpha \log t.}
#' 
#' @param formula A formula containing a \code{\link{Surv}} object, should only
#' contain one categorical predictor, or a set of indicators describing only
#' one predictor.
#' @param data Data set.
#' @param labels A vector containing labels for the plotted lines.
#' @return Produces a plot of log Time vs. log Estimated Cumulative Hazard for
#' each level of the predictor (similarly to what can be obtained using
#' \code{\link{plot.survfit}} and the \code{fun = "cloglog"} option), as well
#' as a data set containing that information.
#' @author Sarah R. Haile
#' @seealso Requires packages \pkg{survival} and \pkg{prodlim}. A similar plot
#' can be produced using \code{\link{plot.survfit}} and the option \code{fun =
#' "cloglog"}.
#' @references Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis:
#' techniques for censored and truncated data}.  Springer.
#' @keywords survival regression
#' @examples
#' 
#' if (requireNamespace("prodlim")) {
#'     data(larynx)
#'     fm <- 
#'     weibullDiag(formula = Surv(time, death) ~ stage, data = larynx,
#'                 labels = c("Stage I", "Stage II","Stage III", "Stage IV"))
#' }
#' 
#' @importFrom prodlim prodlim
#' @export
weibullDiag <- function (formula, data, labels =  NULL) 
{
    ## 'formula' and 'data' checked in prodlim::prodlim()
    np <- prodlim::prodlim(formula = formula, data = data)
    if(is.null(labels))
        labels <- rownames(np$X)
    else
        stopifnot(is.character(labels), length(labels) == nrow(np$X))
    cols <- NA
    for (i in 1:length(np$size.strata)) {
        cols <- c(cols, rep(i, np$size.strata[i]))
    }
    cols <- cols[!is.na(cols)]
    window <- c(-0.5, 0.5)
    y.range <- range(log(-log(np$surv)), finite = TRUE)
    x.range <- log(range(np$time))
    plot(0, 0, type = "n", ylim = y.range + window, xlim = x.range + window,
         xlab = "Log Survival Time", ylab = "Log Cumulative Hazard", 
         main = "Weibull Regression\nDiagnostic Plot")
    for (i in 1:length(np$size.strata)) {
        t <- np$time[cols == i]
        s <- np$surv[cols == i]
        h <- np$haz[cols == i]
        points(log(t), log(-log(s)), col = i, pch = i, lty = i, lwd = 2, type = "b")
    }
    legend("topleft", legend = labels, col = 1:length(np$size.strata), 
           pch = 1:length(np$size.strata), lwd = 2, lty = i:length(np$size.strata))
    pCH <- list(x = log(np$time), y = log(-log(np$surv)), strata = cols, 
                stata.def = np$X)
}
