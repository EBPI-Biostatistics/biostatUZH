#' Confidence intervals for weighted kappa and m >= 2 raters
#'
#' Compute confidence intervals for the coefficient of
#' agreement for two nominal or ordered variables and two or more raters.
#'
#' @param dat data.frame that contains the ratings as columns.
#' @param type Defines the type of confidence interval that is computed. If equal to "Cohen", then Cohen's
#' unweighted \eqn{kappa} is computed, i.e. ratings are assumed to be nominal. If not equal to
#' "Cohen", the weighted version for ordered ratings is computed.
#' @param weights Define weights to be used if ordered ratings are compared. Only used if \code{type != "Cohen"}.
#' @param M Number of bootstrap samples to be generated.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A list containing:
#' \item{n}{Number of observations used to compute confidence intervals.}
#' \item{kappa}{Computed \code{kappa}.}
#' \item{boot.quant}{Confidence interval based on quantiles of the bootstrap distribution.}
#' @references
#' Conger, A.J. (1980), \emph{Integration and generalisation of Kappas for multiple raters},
#' Psychological Bulletin, 88, 322-328.
#' @seealso This function basically implements the example given for
#' \code{lkappa} in package \pkg{psy}.
#' @details
#' This function computes bootstrap confidence intervals for an
#' unweighted or weight \eqn{kappa} coefficient,
#' based on all pairwise complete observations in \code{dat}.
#' @examples
#' if (requireNamespace("psy")) {
#'   ## example comparable to that when called ?lkappa
#'   data("expsy", package = "psy")
#'   set.seed(1)
#'   confIntKappa(dat = expsy[,c(11,13,15)], type = "not Cohen", weights = "absolute",
#'                M = 200, conf.level = 0.95)
#' }
#' @keywords htest
#' @export
confIntKappa <- function(dat, type = "not Cohen",
                         weights = c("absolute", "squared")[1],
                         M = 1000, conf.level = 0.95)
{
## =================================================================
## Fleiss' kappa for m raters according to Conger (1980), allows
## for weighting also in the case of m > 2 raters
## computed using the function lkappa in package 'psy'. 
## For bootstrap confidence interval, package 'boot' is needed.
##
## Input:
##    - dat:     m * n matrix of ratings (m subjects, n raters)
##    - type:    If = "Cohen", then Cohen's unweighted kappa is computed, i.e.
##               ratings are assumed to be nominal. If != "Cohen", the weighted
##               version for ordered ratings is computed.
##    - weights: "absolute" or "squared"
##
## =================================================================

    if (!requireNamespace("psy")) stop("requires psy::lkappa()")
    
    alpha <- 1 - conf.level
    
    ## compute number of complete observations
    n <- sum(apply(is.na(dat), 1, sum) == 0)
    
    ## compute kappa
    k <- psy::lkappa(dat, type = type, weights = weights)
    
    ## generate bootstrap confidence interval
    if ((type == "Cohen") & (weights == "absolute")){kappam.boot <- function(data, x){psy::lkappa(r = data[x, ], type = "Cohen", weights = "absolute")}}
    if ((type == "Cohen") & (weights == "squared")){kappam.boot <- function(data, x){psy::lkappa(r = data[x, ], type = "Cohen", weights = "squared")}}
    if ((type != "Cohen") & (weights == "absolute")){kappam.boot <- function(data, x){psy::lkappa(r = data[x, ], type = "not Cohen", weights = "absolute")}}
    if ((type != "Cohen") & (weights == "squared")){kappam.boot <- function(data, x){psy::lkappa(r = data[x, ], type = "not Cohen", weights = "squared")}}
    res <- boot(data = dat, statistic = kappam.boot, R = M)
    
    ## quantile of bootstrap samples
    quantil <- quantile(res$t, c(alpha / 2, 1 - alpha / 2)) 
    
    ## adjusted bootstrap percentile (BCa) confidence interval (better)
                                        #adj.boot <- boot.ci(res, conf = conf.level, type = "bca")$bca[4:5]    
    
    
                                        # package psy drops levels that are not present in the data for
                                        # ..function "ckappa" & "wkappa", this can lead to wrong calculations
    if(weights == "squared"){
        message((paste0("Caution, used levels in weighted Kappa: ",
                        paste0(levels(as.factor(c(as.character(dat[, 1]),
                                                  as.character(dat[, 2])))), 
                               collapse = ", " ))))
    }
    
    ## generate output
    res <- list("n" = n, "kappa" = k, "boot.quant" = quantil)#, "adj.boot" = adj.boot)
    return(res)
}
