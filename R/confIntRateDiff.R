#' Compute confidence interval for a rate difference
#' 
#' Two methods to compute a confidence interval for a rate difference based on
#' Wald and Wilson confidence intervals for the individual rates are provided.
#' 
#' 
#' @param x vector of length 2, number of events in each group.
#' @param t vector of length 2, total observation time in each group.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries: \item{rd}{Estimated rate difference.}
#' \item{CIs}{Dataframe containing confidence intervals for the rate
#' difference.}
#' @author Leonhard Held
#' @seealso \code{\link{wilsonRate}}
#' @references Held, L., Rufibach, K. and Seifert, B. (2013). Medizinische
#' Statistik - Konzepte, Methoden, Anwendungen.  Section 8.2.
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' t <- c(100, 120)
#' confIntRateDiff(x, t)$CIs
#' 
confIntRateDiff <- function(x, t, conf.level = 0.95){

    ci1 <-  wilsonRate(x[1], t[1], conf.level = conf.level)
    ci2 <-  wilsonRate(x[2], t[2], conf.level = conf.level)

    diff <- matrix(ci1[2] - ci2[2])
    se.diff <- sqrt(x[1]/t[1]^2+x[2]/t[2]^2)
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


