confIntDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95)
{
    stopifnot(is.wholenumber(tp), is.wholenumber(fp),
              is.wholenumber(tn), is.wholenumber(fn),  conf.level<1,
              conf.level>0)

    res <- data.frame(matrix(NA, ncol = 4))
    colnames(res) <- c("type", "lower", "estimate", "upper")

    res[1, 2:4] <- wilson(x=tp, n=tp+fn, conf.level = conf.level)
    res[2, 2:4] <- wilson(x=tn, n=tn+fp, conf.level = conf.level)
    res[3, 2:4] <- confIntRiskRatio(x=c(tp,fp), n=c(tp+fn, fp+tn), conf.level = conf.level)
    res[4, 2:4] <- confIntRiskRatio(x=c(fn,tn), n=c(tp+fn, tn+fp), conf.level = conf.level)

    res[, 1] <- c("Sensitivity", "Specificity", "LRplus", "LRminus")

    return(res)

}
