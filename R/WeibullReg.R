#' Weibull Regression for Survival Data
#' 
#' \code{WeibullReg} performs Weibull regression using the
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
#' in the documentation for \code{\link{ConvertWeibull}}.
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
#' @seealso Requires the package \pkg{survival}. This function depends on
#' \code{\link{ConvertWeibull}}. See also \code{\link{survreg}}.
#' @references Carroll, K. (2003).  On the use and utility of the Weibull model
#' in the analysis of survival data. \emph{Controlled Clinical Trials},
#' \bold{24}, 682--701.
#' 
#' Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis: techniques
#' for censored and truncated data}.  Springer.
#' @keywords survival regression
#' @examples
#' 
#' data("larynx")
#' WR <- WeibullReg(Surv(time, death) ~ factor(stage) + age, data=larynx)
#' WR
#' 
#' @export
WeibullReg <- function (formula, data = parent.frame(), conf.level = 0.95) 
{
    m <- survreg(formula, data, dist = "weibull")
    mle <- ConvertWeibull(m, conf.level)
    return(list(formula = formula, coef = mle$vars, HR = mle$HR, 
        ETR = mle$ETR, summary = summary(m)))
}
