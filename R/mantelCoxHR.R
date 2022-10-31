#' Mantel-Cox estimator of the hazard ratio
#' 
#' Estimate the hazard ratio and compute a confidence interval for it via a
#' special application of the Mantel-Haenszel procedure. A separate 2 x 2 table
#' is constructed for each event time. The underlying assumption is that the
#' hazard ratio is constant over the follow-up period.
#' 
#' 
#' @param time Event times, censored or observed.
#' @param event Censoring indicator, 1 for event, 0 for censored.
#' @param group Factor with two levels, e.g. treatment group.
#' @param conf.level Significance level for confidence interval for hazard
#' ratio.
#' @return \item{mantelCox.hr}{Hazard ratio estimate.}
#' \item{ci.hr}{Wald confidence interval at the level specified by \code{conf.level}.}
#' \item{p.val.logrank}{\eqn{p}-value of logrank test for a comparison of the
#' survival curves between the two groups.}
#' \item{coxph.hr}{Hazard ratio estimated via Cox-regression.}
#' @note In general the Mantel-Cox estimate and the hazard ratio
#' estimate received via Cox-regression do not coincide.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Kirkwood, B.R. and Sterne, J.A.C. (2003). \emph{Essential
#' Medical Statistics}. Blackwell Science. See p. 283 ff.
#' @keywords htest survival
#' @examples
#' 
#' ## load Acute Myelogenous Leukemia survival data
#' library(survival)
#' data("leukemia", package = "survival")
#' 
#' mantelCoxHR(time = leukemia[, 1],
#'             event = leukemia[, 2],
#'             group = as.factor(leukemia[, 3]))
#'
#' @export
mantelCoxHR <- function(time, event, group, conf.level = 0.95){

    ## 'time', 'event', 'group' are checked in 'Surv()' and 'survfit()'
    stopifnot(is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)

    alpha <- 1 - conf.level
    group <- as.factor(group)
    s <- summary(survival::survfit(survival::Surv(time, event) ~ group, conf.type = "none") )
    dat <- data.frame(s$time, s$strata, s$n.risk, s$n.event) 
    dat <- dat[order(s$time), ] 
    lev <- levels(dat$s.strata)
    res <- matrix(NA, ncol = 4, nrow = nrow(dat))
    
    for (i in seq_len(nrow(dat))){
        tmp <- dat[i, ]
        if (tmp[2] == lev[1]){res[i, 1:2] <- c(tmp[[3]], tmp[[4]])}
        if (tmp[2] == lev[2]){res[i, 3:4] <- c(tmp[[3]], tmp[[4]])}
    }
    
    ## res2 is the table on p. 285 in Kirkwood & Sterne
    res2 <- res 
    n2 <- length(res2[, 1])
    for (i in 1:n2){
        g1 <- res2[i, 1:2]
        if (is.na(g1[1]) == TRUE){
            res2[i, 1:2] <- c(sum(time[group == levels(group)[1]] >= dat$s.time[i], na.rm = TRUE), 0)
        }
        
        g2 <- res2[i, 3:4]
        if (is.na(g2[1]) == TRUE){
            res2[i, 3:4] <- c(sum(time[group == levels(group)[2]] >= dat$s.time[i], na.rm = TRUE), 0)
        }
    } 
    dimnames(res2)[[2]] <- c("n0i", "d0i", "n1i", "d1i")
    
    ## compute Mantel-Cox HR
    h0i <- res2[, "n0i"] - res2[, "d0i"] 
    h1i <- res2[, "n1i"] - res2[, "d1i"]
    ni <- res2[, "n0i"] + res2[, "n1i"]
    di <- res2[, "d0i"] + res2[, "d1i"]
    hi <- h0i + h1i
    Q <- sum(res2[, "d1i"] * h0i / ni)
    R <- sum(res2[, "d0i"] * h1i / ni)
    U <- sum(res2[, "d1i"] - di * res2[, "n1i"] / ni)
    V <- sum(di * hi * res2[, "n1i"] * res2[, "n0i"] / (ni ^ 2 * (ni - 1)), na.rm = TRUE)
    hr <- Q / R
    
    ## compute confidence interval
    se.log.hr <- sqrt(V / (Q * R))
    ci.hr <- exp(log(hr) + stats::qnorm(c(alpha / 2, 1 - alpha / 2)) * se.log.hr)
    
    ## p-value logrank test
    chi2 <- U ^ 2 / V
    p.val <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
    
    ## compare to HR computed from Cox-regression 
    fit2 <- survival::coxph(survival::Surv(time, event) ~ group)
    cfit2 <- as.numeric(stats::coef(fit2))
    hr.cox <- exp(cfit2)
    
    list("mantelCox.hr" = hr, "ci.hr" = ci.hr, "p.val.logrank" = p.val, "coxph.hr" = hr.cox)
}
