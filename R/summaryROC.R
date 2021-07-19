#' Function to compute ROC curve and an asymptotic confidence interval for AUC
#' 
#' This function computes for values of a continuous variable of a group of
#' cases and a group of controls the ROC curve. Additionally, the AUC including
#' an asymptotic confidence interval is provided.
#' 
#' 
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param conf.level Confidence level of confidence interval.
#' @return \item{x.val}{1-specificity of the test, so the values on the
#' \eqn{x}-axis of a ROC plot.} \item{y.val}{Sensitivity of the test, so the
#' values on the \eqn{y}-axis of a ROC plot.} \item{ppvs}{Positive predictive
#' values for each cutoff.} \item{npvs}{Negative predictive values for each
#' cutoff.} \item{cutoffs}{Cutoffs used (basically the pooled marker values of
#' cases and controls).} \item{res.mat}{Collects the above quantities in a
#' matrix, including Wilson confidence intervals, computed at at confidence
#' level \code{conf.level}.} \item{auc}{Area under the ROC curve. This is equal
#' to the value of the Mann-Whitney test statistic.} \item{auc.var}{Variance of
#' AUC.} \item{auc.var.norm}{Variance of AUC if data is assumed to come from a
#' bivariate normal distribution.} \item{lowCI}{Lower limit of Wald confidence
#' interval.} \item{upCI}{Upper limit of Wald confidence interval.}
#' \item{logitLowCI}{Lower limit of a Wald confidence interval received on
#' logit scale.} \item{logitUpCI}{Upper limit of a Wald confidence interval
#' received on logit scale.}
#' @note The confidence intervals are only valid if observations are
#' \emph{independent}.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}. Part of the
#' function was derived from code by Andrea Riebler.
#' @seealso Similar functionality is provided in the package \pkg{ROCR}.
#' However, this latter package offers no computation of confidence intervals.
#' @references The original reference for computation of confidence intervals
#' is:
#' 
#' Hanley, J.A. and McNeil, B.J. (1982). The meaning and use of the area under
#' the curve. \emph{Radiology}, \bold{143}, 29--36.
#' 
#' See also:
#' 
#' Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for
#' classification and prediction}. Oxford: Oxford University Press.
#' @keywords htest
#' @examples
#' 
#' 
#' set.seed(1977)
#' ns <- c(50, 40)
#' truth <- c(rep(0, ns[1]), rep(1, ns[2]))
#' estimates <- c(rnorm(ns[1]), rnorm(ns[2], mean = 0.5, sd = 1.5))
#' cases <- estimates[truth == 1]
#' controls <- estimates[truth == 0]
#' res <- summaryROC(cases, controls, conf.level = 0.95)
#' 
#' # display results
#' res
#' res$res.mat
#' 
#' # plot ROC curve
#' plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
#'     xlab = "1 - specificity", ylab = "sensitivity", pty = 's')
#' segments(0, 0, 1, 1, lty = 2)
#' lines(res$x.val, res$y.val, type = 'l', col = 2, lwd = 2, lty = 2)
#' 
summaryROC <- function(cases, controls, conf.level = 0.95){
    
    # Compute:
    #
    #   - ROC curve
    #   - AUC, including confidence interval
    #
    # Input
    #   - cases: Marker-values of cases
    #   - controls: Marker-values of controls
    #
    # Inspired by code of A. Riebler. CI computed according to Pepe (2003), p. 105 ff.
    #
    # Kaspar Rufibach, September 2008
    
    alpha <- 1 - conf.level

    n1 <- length(controls)
    n2 <- length(cases)
    cutoffs <- c(Inf, rev(sort(unique(c(cases, controls)))))
    n <- length(cutoffs)
    tps <- rep(NA, n)
    tns <- tps; fps <- tps; fns <- tps; n1tmp <- tps; n3tmp <- tps

    for (i in 1:n){
        cut <- cutoffs[i]
        
        fps[i] <- sum(controls >= cut)    
        tps[i] <- sum(cases >= cut)
        tns[i] <- n1 - fps[i]    
        fns[i] <- n2 - tps[i]   
        
        # auxiliary quantities to compute varAUC
        n1tmp[i] <- sum(controls == cut)
        n3tmp[i] <- sum(cases == cut)
    }

    x.val <- fps / (tns + fps)     #1 - tns / (tns + fps)
    y.val <- tps / (tps + fns)
    
    ppvs <- tps / (tps + fps)
    npvs <- tns / (fns + tns)
    
    # compute confidence intervals
    spec.ci <- matrix(NA, nrow = n, ncol = 2)
    sens.ci <- spec.ci
    ppvs.ci <- spec.ci
    npvs.ci <- spec.ci
    
    for (i in 1:n){
        spec.ci[i, ] <- wilson(x = tns[i], n = (tns + fps)[i], conf.level = conf.level)[c(1, 3)]
        sens.ci[i, ] <- wilson(x = tps[i], n = (tps + fns)[i], conf.level = conf.level)[c(1, 3)]
        ppvs.ci[i, ] <- wilson(x = tps[i], n = (tps + fps)[i], conf.level = conf.level)[c(1, 3)]
        npvs.ci[i, ] <- wilson(x = tns[i], n = (fns + tns)[i], conf.level = conf.level)[c(1, 3)]
    }
    
    # compute q2 and q1
    q2 <- sum(n3tmp * (tns ^ 2 + tns * n1tmp + 1/3 * n1tmp ^ 2)) / (n2 * n1 ^ 2)
    q1 <- sum(n1tmp *(tps ^ 2 + tps * n3tmp + 1/3 * n3tmp ^ 2)) / (n1 * n2 ^ 2)

    # estimate auc according to Hanley and McNeil (1982)
    auc <- sum(n1tmp * tps + 0.5 * n1tmp * n3tmp) / (n1 * n2)

    # estimate AUC as test statistic of Wilcoxon test
    auc <- as.numeric(wilcox.test(cases, controls, exact = FALSE)$statistic / (n1 * n2))
    auc.var <- (n1 * n2) ^ (-1) * (auc * (1 - auc) + (n2 - 1) * (q1 - auc ^ 2) + (n1 - 1) * (q2 - auc ^ 2))

    # compute variance of AUC when assuming a normal model for measurements
    q1norm <- auc / (2 - auc)
    q2norm <- 2 * auc ^ 2 / (1 + auc)
    auc.var.norm <- (n1 * n2) ^ (-1) * (auc * (1 - auc) + (n2 - 1) * (q1norm - auc ^ 2) + (n1 - 1) * (q2norm - auc ^ 2))

    # compute confidence intervals
    # on original scale
    auc.se <- sqrt(auc.var)
    q <- qnorm(1 - alpha / 2)
    lowCI <- auc - q * auc.se
    upCI <- auc + q * auc.se
    
    logitAuc <- log(auc / (1 - auc))

    # using the Delta rule: d / d auc logit(auc))* auc.se
    logitAucSE <- auc.se / (auc * (1 - auc))
    logitLowCI <- logitAuc - q * logitAucSE
    logitUpCI <- logitAuc + q * logitAucSE

    # backtransformation
    logLowCI <- 1 / (1 + exp(- logitLowCI))
    logUpCI <- 1 / (1 + exp(- logitUpCI))
    
    # matrix containing most important stuff, incl. cutoff and confidence intervals
    res.mat <- cbind(cutoffs, 1 - x.val, spec.ci, y.val, sens.ci, ppvs, ppvs.ci, npvs, npvs.ci)
    colnames(res.mat) <- c("cutoff", "specificity", "CIspeclow", "CIspecup", "sensitivity", "CIsenslow", "CIsensup", "PPV", "CIPPVlow", "CIPPVup", "NPV", "CINPVlow", "CINPVup")

    res <- list("x.val" = x.val, "y.val" = y.val, "ppvs" = ppvs, "npvs" = npvs, "cutoffs" = cutoffs, "res.mat" = res.mat, "auc" = auc, 
        "auc.var" = auc.var, "auc.var.norm" = auc.var.norm, "lowCI" = lowCI, "upCI" = upCI, "logitLowCI" = logLowCI, "logitUpCI" = logUpCI)
        
    return(res)
}
