confIntAUC <- function(cases, controls, conf.level = 0.95){
    
    # Compute:
    #   - AUC, including confidence interval, including one on logit scale
    # Input
    #   - cases: Marker-values of cases
    #   - controls: Marker-values of controls
    #                                    #
    # Using parts of a previous implementation by Kaspar Rufibach
                                        # Leonhard Held, September 2016
    
    ##     # we follow the second method described in Altman et al, p. 113, 
    ## standardErrorAUC <- function(cases, controls){
    ##     ncases <- length(cases)
    ##     ncontrols <- length(controls)
    ##     # non-disease placement values of cases and controls
    ##     C <- rep(NA, ncases)
    ##     R <- rep(NA, ncontrols)
        
    ##     for(i in 1:ncases)
    ##         C[i] <- mean(as.numeric(controls<cases[i])+0.5*as.numeric(controls==cases[i]))
    ##     for(j in 1:ncontrols)
    ##         R[j] <- mean(as.numeric(cases>controls[j])+0.5*as.numeric(cases==controls[j]))
    ##     auc.se <- sqrt((var(R)/ncontrols + var(C)/ncases))
    ##     return(auc.se)
    ## }
    
    # estimate AUC as normalized test statistic of Wilcoxon test
    ncontrols <- length(controls)
    ncases <- length(cases)
    auc <- as.numeric(wilcox.test(cases, controls, exact = FALSE)$statistic / (ncases * ncontrols))
    auc.se <- standardErrorAUC(cases, controls)

    # compute confidence intervals
    # on original scale
    z <- qnorm((1 + conf.level) / 2)
    lower <- auc - z * auc.se
    upper <- auc + z * auc.se

    # on logit scale
    logitAuc <- log(auc / (1 - auc))
    logitAucSE <- auc.se / (auc * (1 - auc))
    logitLowCI <- logitAuc - z * logitAucSE
    logitUpCI <- logitAuc + z * logitAucSE

    # backtransformation
    lowerLogit <- 1 / (1 + exp(- logitLowCI))
    upperLogit <- 1 / (1 + exp(- logitUpCI))
    
    res <- data.frame(matrix(NA, ncol = 4))
    colnames(res) <- c("type", "lower", "AUC", "upper")
    res[1, 2:4] <- c(lower, auc, upper)
    res[2, 2:4] <- c(lowerLogit, auc, upperLogit)
    res[, 1] <- c("Wald", "logit Wald")
    return(res)
}
