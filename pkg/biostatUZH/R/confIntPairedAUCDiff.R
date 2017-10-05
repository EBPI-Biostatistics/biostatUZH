
## computes confidence interval for the difference in AUC of two paired tests
## cases is a m x 2 matrix or data frame
## controls is a n x 2 matrix or data frame
## each row corresponds to measurements from one case/control

confIntPairedAUCDiff <- function(cases, controls, conf.level = 0.95){
    stopifnot(is.matrix(cases) | is.data.frame(cases), is.matrix(controls) | is.data.frame(controls))

    # estimate AUC as normalized test statistic of Wilcoxon test
    ncontrols <- nrow(controls)
    ncases <- nrow(cases)
    auc <- numeric()
    for(k in 1:2)
        auc[k] <- as.numeric(wilcox.test(cases[,k], controls[,k], exact = FALSE)$statistic / (ncases * ncontrols))
    auc.diff <- auc[1]-auc[2]
    auc.diff.se <- standardErrorAUCDiff(cases, controls)

    # compute confidence intervals
    # on original scale
    z <- qnorm((1 + conf.level) / 2)
    lower <- auc.diff - z * auc.diff.se
    upper <- auc.diff + z * auc.diff.se

    res <- data.frame(matrix(NA, ncol = 4))
    colnames(res) <- c("outcome", "lower", "estimate", "upper")
    res[1, 2:4] <- confIntAUC(cases[,1], controls[,1])[2,2:4] ## avoids overshoot
    res[2, 2:4] <- confIntAUC(cases[,2], controls[,2])[2,2:4] ## avoids overshoot
    res[3, 2:4] <- c(lower, auc.diff, upper)
    res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")

    return(res)
}
