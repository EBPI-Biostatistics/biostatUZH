#' Compute confidence interval for a Poisson rate via different methods
#' 
#' Compute a confidence interval for a Poisson rate using several methods. The
#' individual methods are also available as separate functions \code{waldRate}
#' and \code{wilsonRate}.
#' 
#' 
#' @aliases confIntRate waldRate wilsonRate
#' @param x Number of events.
#' @param t Total observation time.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries: \item{p}{Estimated rate.}
#' \item{CIs}{Dataframe containing the estimated confidence intervals.}
#' @author Leonhard Held
#' @seealso Functions for some of the intervals provided here are available in
#' \pkg{Hmisc} (see the examples).
#' @references Held, L., Rufibach, K. and Seifert, B. (2013). Medizinische
#' Statistik - Konzepte, Methoden, Anwendungen.  Section 8.2.
#' @examples
#' 
#' ## Calculate confidence bounds for a Poisson rate by different methods. 
#' x <- 1
#' t <- 3
#' ci <- confIntRate(x, t)$CIs
#' ci
#' p <- confIntRate(x, t)$rate
#' 
#' plot(0, 0, type = 'n', ylim = c(0, 3), xlim = c(0, 3), xlab = 'p', 
#'     ylab = '', yaxt = 'n')
#' lines(ci[1, 2:3], c(1, 1))
#' lines(ci[2, 2:3], c(2, 2))
#' points(p, 1, pch=19)
#' points(p, 2, pch=19)
#' text(0.5, 0.85, 'wald')
#' text(0.5, 1.85, 'wilson')
#' 
confIntRate <- function(x, t, conf.level = 0.95)
{
    stopifnot(is.wholenumber(x), (t>0),  conf.level<1, conf.level>0)

    res <- data.frame(matrix(NA, ncol = 3))
    colnames(res) <- c("type", "lower", "upper")

    res[1, 2:3] <- waldRate(x, t, conf.level = conf.level)[c(1, 3)]
    res[2, 2:3] <- wilsonRate(x, t, conf.level = conf.level)[c(1, 3)]

    res[, 1] <- c("Wald", "Wilson")

    res <- list("rate" = x / t, "CIs" = res)
    return(res)

}
