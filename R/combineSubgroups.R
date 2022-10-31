#' Combines summary data across subgroups
#' 
#' Combines mean estimates, standard errors, and sample sizes of
#' several treatment and placebo groups into an overall treatment effect
#' and its standard error. Assumes normally distributed estimates. 
#' 
#' @param n Vector with sample sizes in each subgroup.
#' @param means Vector with sample means in each subgroup.
#' Has to be of same length as \code{n}.
#' @param se Vector with sample standard errors in each subgroup.
#' Has to be of same length as \code{n}.
#' @param treatment Index with intagers between 1 and \code{length(n)}.
#' Indicates which entries of each vector belong to the
#' treatment group.
#' @return Overall difference in means with a standard error.
#' @author Leonhard Held
#' @examples
#' 
#' combineSubgroups(n = c(10, 20, 30, 40), means = c(12, 11, 10, 9),
#'                  se = c(3, 4, 3, 4), treatment = c(2, 4))
#' 
#' @export
combineSubgroups <- function(n, means, se, treatment){

    stopifnot(is.numeric(n), length(n) >= 2,
              is.finite(n), is.wholenumber(n),
              is.numeric(means), length(n) == length(means),
              is.finite(n),
              is.numeric(se), length(n) == length(se),
              is.finite(se),
              is.numeric(treatment), length(treatment) <= n,
              is.finite(treatment), is.wholenumber(treatment),
              1 <= treatment, treatment <= length(n))
    
    ## sample sizes for intervention and placebo
    n.trtm <- n[treatment]
    n.plac <- n[-treatment]
    ## mean response in both groups
    theta.trtm <- stats::weighted.mean(means[treatment], w = n.trtm)
    theta.plac <- stats::weighted.mean(means[-treatment], w = n.plac)
    ## overall treatment effect
    theta <- theta.trtm - theta.plac
    ## within-group variance
    varw.trtm <- stats::weighted.mean(se[treatment]^2, w = n.trtm)
    varw.plac <- stats::weighted.mean(se[-treatment]^2, w = n.plac)
    ## between-group variance
    varb.trtm <- stats::weighted.mean((means[treatment] - theta.trtm)^2, w = n.trtm)
    varb.plac <- stats::weighted.mean((means[-treatment] - theta.plac)^2, w = n.plac)
    ## total variance
    var.trtm <- varw.trtm + varb.trtm
    var.plac <- varw.plac + varb.plac
    ## standard error of overall treatment effect
    se.theta <- sqrt(var.trtm / sum(n.trtm) + var.plac / sum(n.plac))
    res <- c(theta, se.theta)
    names(res) <- c("effect estimate", "standard error")
    res
}
 
