#' Confidence interval for a given effect
#' 
#' Prints the Wald confidence interval and the p-value for a given effect and
#' standard error.
#' 
#' 
#' @param theta the effect size.
#' @param se.theta the standard error.
#' @param conf.level Confidence level. Default is 95\%.
#' @param FUN Transformation to applied to effect size and confidence interval.
#' Default is the identity function.
#' @param digits number of digits for the effect.
#' @return Prints effect and associated Wald confidence interval and p-value.
#' @author Leonhard Held
#' @examples
#' 
#' printResults(1.7, 2.1)
#' printResults(1.7, 2.1, FUN=exp)
#' 
printResults <- function(theta, se.theta, conf.level=0.95, FUN=identity, digits=3){
    z <- qnorm( (1 + conf.level) / 2)
    ci.l <- theta - z * se.theta
    ci.u <- theta + z * se.theta
    p <-  2 * pnorm(abs(theta) / se.theta, lower.tail=FALSE)
    effect <- round(FUN(theta), digits=digits)
    ci <- formatCI(c(FUN(ci.l), FUN(ci.u)), digits=digits, text="english")
    p <- formatPval(p)
    res <- c(effect, ci, p)
    names(res) <- list("Effect", 
                       sprintf('%d%% Confidence Interval', conf.level * 100),
                       "P-value")
    print(t(res), quote=FALSE)
}
