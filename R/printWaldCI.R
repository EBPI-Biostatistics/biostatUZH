#' Confidence Wald confidence interval for a given effect
#' 
#' Takes an effect and its standard error to calculate the Wald
#' confidence interval and the p-value.
#' 
#' @param theta the effect size.
#' @param se.theta the standard error.
#' @param conf.level Confidence level. Default is 95\%.
#' @param FUN Transformation to applied to effect size and confidence interval.
#' Default is the identity function.
#' @param digits number of digits for the effect.
#' @return Prints and invisiblie returns the effect, the Wald confidence interval, and the p-value.
#' @author Leonhard Held
#' @examples
#' 
#' printWaldCI(theta = 1.7, se.theta = 2.1)
#' printWaldCI(theta = 1.7, se.theta = 2.1, FUN = exp)
#' 
#' @export
printWaldCI <- function(theta, se.theta, conf.level = 0.95, FUN = identity, digits = 3){
    stopifnot(is.numeric(theta), length(theta) == 1,
              is.finite(theta),
              is.numeric(se.theta), length(se.theta) == 1,
              is.finite(se.theta), se.theta > 0,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1,
              is.function(FUN),
              is.numeric(digits), length(digits) == 1,
              is.finite(digits), is.wholenumber(digits), 0 <= digits)

    z <- qnorm((1 + conf.level) / 2)
    ci.l <- theta - z * se.theta
    ci.u <- theta + z * se.theta
    p <-  2 * pnorm(abs(theta) / se.theta, lower.tail = FALSE)
    effect <- round(FUN(theta), digits = digits)
    ci <- formatCI(c(FUN(ci.l), FUN(ci.u)), digits = digits, text = "english")
    p <- formatPval(p)
    res <- c(effect, ci, p)
    names(res) <- list("Effect", 
                       sprintf('%d%% Confidence Interval', conf.level * 100),
                       "P-value")
    res <- t(res)
    print(x = res, quote = FALSE)
    invisible(res)
}
