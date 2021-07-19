#' Confidence interval for a risk differences,
#' 
#' Two methods to compute a confidence interval for a risk difference based on
#' Wald and Wilson confidence intervals for the individual risks are provided.
#' 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A list with the entries: \itemize{ \itemrdEstimated risk difference.
#' \itemCIsDataframe containing confidence intervals for the risk difference. }
#' @author Leonhard Held
#' @seealso \code{\link{wilson}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntOddsRatio}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 6
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntRiskDiff(x, n)$CIs
#' 
confIntRiskDiff <- function(x, n, conf.level = 0.95){

    ci1 <-  wilson(x[1], n[1], conf.level = conf.level)
    ci2 <-  wilson(x[2], n[2], conf.level = conf.level)

    diff <- matrix(ci1[2] - ci2[2])
    se.diff <- sqrt((ci1[2] * (1 - ci1[2])) / n[1] + (ci2[2] * (1 - ci2[2])) / n[2])  
    z <- qnorm((1 + conf.level) / 2)
    wald.lower <- diff - z * se.diff
    wald.upper <- diff + z * se.diff
    score.lower <- diff - sqrt((ci1[2] - ci1[1])^2 + (ci2[3] - ci2[2])^2)
    score.upper <- diff + sqrt((ci2[2] - ci2[1])^2 + (ci1[3] - ci1[2])^2)

    result <- matrix(ncol=2, nrow=2)
    result[,1] <- c(wald.lower, score.lower)
    result[,2] <- c(wald.upper, score.upper)

    out <- data.frame(type = c("Wald", "Wilson"), result)
    names(out) <- c("type", "lower", "upper")

    ret.list <- list("rd" = diff, "CIs" = out)
    return(ret.list)

}




#' Confidence interval for a risk ratio
#' 
#' Provides a confidence interval for a risk ratio. The method is based on a
#' Wald interval for the log risk ratio. Used by
#' \code{\link{confIntDiagnostic}}.
#' 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A vector containing the risk ratio and the limits of the confidence
#' interval.
#' @author Leonhard Held
#' @seealso \code{\link{confIntDiagnostic}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntOddsRatio}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 7
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntRiskRatio(x, n)
#' 
confIntRiskRatio <- function(x, n, conf.level = 0.95){
    stopifnot(length(x)==2, length(n)==2, is.wholenumber(x), is.wholenumber(n), (x>0), (x<n),
              conf.level<1, conf.level>0)
    Risk <- x / n
    RiskRatio <- Risk[1] / Risk[2]
    se.log.RiskRatio <- sqrt(sum(1 / x) - sum(1 / n))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * se.log.RiskRatio)
    wald.lower <- RiskRatio / EF
    wald.upper <- RiskRatio * EF

    res <- c("lower"=wald.lower, "Risk Ratio"=RiskRatio, "upper"=wald.upper)
    return(res)

}




#' Confidence interval for an odds ratio
#' 
#' Provides a confidence interval for an odds ratio. The method is based on a
#' Wald interval for the log odds ratio. Used by
#' \code{\link{confIntDiagnostic}}.
#' 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A vector containing the odds ratio and the limits of the confidence
#' interval.
#' @author Leonhard Held
#' @seealso \code{\link{confIntDiagnostic}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntRiskDiff}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 7
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntOddsRatio(x, n)
#' 
confIntOddsRatio <- function(x, n, conf.level = 0.95){
    stopifnot(length(x)==2, length(n)==2, is.wholenumber(x), is.wholenumber(n), (x>0), (x<n),
              conf.level<1, conf.level>0)
    y <- n - x
    Odds <- x / y
    OddsRatio <- Odds[1] / Odds[2]
    se.log.OddsRatio <- sqrt(sum(1 / x) + sum(1 / y))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * se.log.OddsRatio)
    wald.lower <- OddsRatio / EF
    wald.upper <- OddsRatio * EF

    res <- c("lower"=wald.lower, "Odds Ratio"=OddsRatio, "upper"=wald.upper)
    return(res)
}


