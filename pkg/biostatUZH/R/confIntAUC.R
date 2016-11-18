confIntAUC <- function(cases, controls, conf.level = 0.95){
    
    # Compute:
    #   - AUC, including confidence interval, including one on logit scale
    # Input
    #   - cases: Marker-values of cases
    #   - controls: Marker-values of controls
    # 
    # Using parts of a previous implementation by Kaspar Rufibach
                                        # Leonhard Held, September 2016
    
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
