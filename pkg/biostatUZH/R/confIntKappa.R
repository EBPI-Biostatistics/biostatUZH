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
#adj.boot <- boot.ci(res, type = "bca")$bca[4:5]        
adj.boot <- NA   

## generate output
res <- list("n" = n, "kappa" = k, "boot.quant" = quantil, "adj.boot" = adj.boot)
return(res)
}







#
