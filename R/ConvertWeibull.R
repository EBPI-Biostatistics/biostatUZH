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
#' function is used by \code{\link{WeibullReg}}.
#' @references Carroll, K. (2003).  On the use and utility of the Weibull model
#' in the analysis of survival data. \emph{Controlled Clinical Trials},
#' \bold{24}, 682--701.
#' 
#' Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis: techniques
#' for censored and truncated data}.  Springer.
#' @keywords survival regression
#' @examples
#' 
#' data(larynx)
#' ConvertWeibull(survreg(Surv(time, death) ~ stage + age, larynx), conf.level = 0.95)
#' 
#' @export
ConvertWeibull <- function (model, conf.level = 0.95)
{
    alpha <- 1 - conf.level
    qa <- qnorm(1 - alpha/2)
    Int.Only <- (nrow(summary(model)$table) == 2)
    sigma <- summary(model)$scale
    mu <- summary(model)$coef[1]
    k <- length(summary(model)$coef) - 1
    if (!Int.Only) {
        gamma <- summary(model)$coef[2:(k + 1)]
    }
    lambda <- exp(-mu/sigma)
    alpha <- 1/sigma
    tmp <- c(lambda, alpha)
    names(tmp) <- c("lambda", "alpha")
    if (!Int.Only) {
        beta <- -gamma/sigma
        tmp <- c(lambda, alpha, beta)
        names(tmp) <- c("lambda", "alpha", names(summary(model)$coef[2:(k + 
            1)]))
    }
    var1 <- summary(model)$var
    var.mu <- diag(var1)[1]
    var.sigma <- var1[(k + 2), (k + 2)] * exp(2 * log(sigma))
    if (!Int.Only) {
        var.gamma <- var1[2:(k + 1), 2:(k + 1)]
        if(k>1) var.gamma <- diag(var.gamma)
        se.gamma <- sqrt(var.gamma)
    }
    cov.mu.sigma <- var1[(k + 2), 1] * sigma
    var.alpha <- var.sigma/(sigma^4)
    var.lambda <- exp(-2 * mu/sigma) * ((var.mu/(sigma^2)) - 
        ((2 * mu/(sigma^3)) * cov.mu.sigma) + (((mu^2)/(sigma^4)) * 
        var.sigma))
    var <- c(sqrt(var.lambda), sqrt(var.alpha))
    if (!Int.Only) {
        cov.gamma.sigma <- var1[2:(k + 1), (k + 2)] * sigma
        var.beta <- (1/(sigma^2)) * (var.gamma - (2 * gamma/sigma) * 
            (cov.gamma.sigma) + (((gamma/sigma)^2) * var.sigma))
        se.beta <- sqrt(var.beta);
        var <- c(sqrt(var.lambda), sqrt(var.alpha), se.beta)
        HR <- cbind(HR = exp(beta), LB = exp(beta - qa * se.beta), 
            UB = exp(beta + qa * se.beta))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
        ETR <- cbind(ETR = exp(gamma), LB = exp(gamma - qa * 
            se.gamma), UB = exp(gamma + qa * se.gamma))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
    }
    tmp1 <- rbind(tmp, var)
    rownames(tmp1) <- c("Estimate", "SE")
    ret <- list(vars = t(tmp1))
    if (!Int.Only) {
        ret$HR = HR
        ret$ETR = ETR
    }
    return(ret)
}


