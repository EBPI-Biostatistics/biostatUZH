#' Combines summary data across subgroups
#' 
#' Combines means and variances of a continuous outcome given in subgroups.
#' 
#' 
#' @param n Vector with sample sizes in each subgroup.
#' @param means Vector with sample means in each subgroup.
#' @param variances Vector with sample variances in each subgroup.
#' @param treatment Indicates which entries of each vector belong to the
#' treatment group.
#' @return Overall difference in means with a standard error.
#' @author Leonhard Held
#' @examples
#' 
#' combineSubgroups(n = c(10, 20, 30, 40), means = c(12, 11, 10, 9),
#'                  variances = c(3, 4, 3, 4), treatment = c(2, 4))
#' 
combineSubgroups <- function(n, means, variances, treatment){

    stopifnot(is.numeric(n), length(n) >= 2,
              is.finite(n), is.wholenumber(n),
              is.numeric(means), length(n) == length(means),
              is.finite(n),
              is.numeric(variances), length(n) == length(variances),
              is.finite(variances),
              is.numeric(treatment), length(n) >= 1,
              is.finite(variances), is.wholenumber(treatment),
              1 <= treatment, treatment <= length(n))
    
    ## sample sizes for intervention and placebo
    n.trtm <- n[treatment]
    n.plac <- n[-treatment]
    ## mean response in both groups
    theta.trtm <- weighted.mean(means[treatment], w = n.trtm)
    theta.plac <- weighted.mean(means[-treatment], w = n.plac)
    ## overall treatment effect
    theta <- theta.trtm - theta.plac
    ## within-group variance
    varw.trtm <- weighted.mean(variances[treatment], w = n.trtm)
    varw.plac <- weighted.mean(variances[-treatment], w = n.plac)
    ## between-group variance
    varb.trtm <- weighted.mean((means[treatment] - theta.trtm)^2, w = n.trtm)
    varb.plac <- weighted.mean((means[-treatment] - theta.plac)^2, w = n.plac)
    ## total variance
    var.trtm <- varw.trtm + varb.trtm
    var.plac <- varw.plac + varb.plac
    ## standard error of overall treatment effect
    se.theta <- sqrt(var.trtm / sum(n.trtm) + var.plac / sum(n.plac))
    res <- c(theta, se.theta)
    names(res) <- c("effect estimate", "standard error")
    res
}
 
