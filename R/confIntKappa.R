#' Confidence intervals for weighted kappa and m >= 2 raters
#' 
#' Compute confidence intervals for the coefficient of agreement for two
#' nominal or ordered variables and two or more raters.
#' 
#' This function computes bootstrap confidence intervals for an unweighted or
#' weight \eqn{kappa} coefficient, based on all pairwise complete observations
#' in \code{data}.
#' 
#' @param data m x n matrix or data.frame containing data from m subjects and n raters.
#' @param type String defining the type of confidence interval.
#' If not equal to "not Cohen" (default), the weighted
#' version for ordered ratings is computed.
#' If "Cohen", Cohen's unweighted \eqn{kappa} is computed, i.e.
#' ratings are assumed to be nominal. 
#' @param weights Define weights to be used if ordered ratings are compared.
#' Can be "squared" (default) or "absolute".
#' Only used if \code{type != "Cohen"}.
#' @param m Number of bootstrap samples to be generated.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A list containing: \item{n}{Number of observations used to compute
#' confidence intervals.} \item{kappa}{Computed \code{kappa}.}
#' \item{boot.quant}{Confidence interval based on quantiles of the bootstrap
#' distribution.}
#' @details Fleiss' kappa for m raters according to Conger (1980), allows
#' for weighting also in the case of m > 2 raters computed using the function
#' \code{\link[psy]{lkappa}} in package \pkg{psy}. 
#' For bootstrap confidence interval the package \pkg{boot} is used.
#' @seealso \code{\link[psy]{lkappa}}, \code{\link[boot]{boot}} 
#' @references Conger, A.J. (1980), \emph{Integration and generalisation of
#' Kappas for multiple raters}, Psychological Bulletin, 88, 322-328.
#' @keywords htest
#' @examples
#' 
#' ## example is similar to example in ?lkappa
#' data("expsy", package = "psy")
#' set.seed(14)
#' confIntKappa(data = expsy[, c(11,13,15)], type = "not Cohen", weights = "absolute",
#'              m = 200, conf.level = 0.95)
#' 
#' @importFrom psy lkappa
#' @importFrom boot boot
#' @export
confIntKappa <- function(data,
                         type = c("not Cohen", "Cohen"),
                         weights = c("squared", "absolute"),
                         m = 1000,
                         conf.level = 0.95)
{
    stopifnot(is.data.frame(data) || is.matrix(data),
              length(data) >= 1, 
              !is.null(type))
    type <- match.arg(type)
    stopifnot(!is.null(weights))
    weights <- match.arg(weights)
    stopifnot(is.numeric(m),
              length(m) == 1,
              is.finite(m),
              is.wholenumber(m),
              1 <= m,
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    alpha <- 1 - conf.level
    
    ## compute number of complete observations
    n <- sum(apply(is.na(data), 1, sum) == 0)
    
    ## compute kappa
    k <- psy::lkappa(data, type = type, weights = weights)

 
    ## generate bootstrap confidence interval
    kappam.boot <- function(data, x) {
        psy::lkappa(r = data[x,], type = type, weights = weights)
    }
    res <- boot(data = data, statistic = kappam.boot, R = m)
    quantil <- quantile(x = res$t, probs = c(alpha / 2, 1 - alpha / 2)) 
    
    if(weights == "squared"){
        message(paste0("Caution, used levels in weighted Kappa: ",
                       paste0(levels(as.factor(c(as.character(data[, 1]),
                                                  as.character(data[, 2])))), 
                               collapse = ", " )))
    }
    
    list("n" = n, "kappa" = k, "boot.quant" = quantil)
}
