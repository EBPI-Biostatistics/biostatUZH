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
