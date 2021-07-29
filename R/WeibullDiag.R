#' Diagnostic Plot of Adequacy of Weibull Distribution
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
#'     fm <- Surv(time, death) ~ stage
#'     WeibullDiag(fm, larynx, labels=c("Stage I", "Stage II","Stage III", "Stage IV"))
#' }
#' 
#' @export
WeibullDiag <- function (formula, data = parent.frame(), labels = rownames(np$X)) 
{
    if (!requireNamespace("prodlim")) stop("requires prodlim::prodlim()")
    np <- prodlim::prodlim(formula, data)
    cols <- NA
    for (i in 1:length(np$size.strata)) {
        cols <- c(cols, rep(i, np$size.strata[i]))
    }
    cols <- cols[!is.na(cols)]
    window <- c(-0.5, 0.5)
    y.range <- range(log(-log(np$surv)), finite = TRUE)
    x.range <- log(range(np$time))
    plot(0, 0, type = "n", ylim = y.range + window, xlim = x.range + 
        window, xlab = "Log Survival Time", ylab = "Log Cumulative Hazard", 
        main = "Weibull Regression\nDiagnostic Plot")
    for (i in 1:length(np$size.strata)) {
        t <- np$time[cols == i]
        s <- np$surv[cols == i]
        h <- np$haz[cols == i]
        points(log(t), log(-log(s)), col = i, pch = i, lty = i, 
            lwd = 2, type = "b")
    }
    legend("topleft", legend = labels, col = 1:length(np$size.strata), 
        pch = 1:length(np$size.strata), lwd = 2, lty = i:length(np$size.strata))
    pCH <- list(x = log(np$time), y = log(-log(np$surv)), strata = cols, 
        stata.def = np$X)
}
