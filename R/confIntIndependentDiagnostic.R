confIntIndependentDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, adjust=FALSE)
{
    stopifnot(is.wholenumber(tp), is.wholenumber(fp),
              is.wholenumber(tn), is.wholenumber(fn),  conf.level<1,
              conf.level>0)
    stopifnot(length(tp)==2, length(fp)==2, length(tn)==2, length(fn)==2)

    se.log.TPF <- sqrt(sum(1/tp)-sum(1/(tp+fn)))
    se.log.FPF <- sqrt(sum(1/fp)-sum(1/(fp+tn)))
    se.log.TNF <- sqrt(sum(1/tn)-sum(1/(tn+fp)))
    se.log.FNF <- sqrt(sum(1/fn)-sum(1/(fn+tp)))

    se.log.rTPF <- sqrt(sum(se.log.TPF^2))
    se.log.rTNF <- sqrt(sum(se.log.TNF^2))
    se.log.rFNF <- sqrt(sum(se.log.FNF^2))
    se.log.rFPF <- sqrt(sum(se.log.FPF^2))
    se.log.rLRp <- sqrt(sum(se.log.TPF^2+se.log.FPF^2))
    se.log.rLRm <- sqrt(sum(se.log.TNF^2+se.log.FNF^2))

    resultA <- confIntDiagnostic(tp[1], fp[1], tn[1], fn[1])
    resultB <- confIntDiagnostic(tp[2], fp[2], tn[2], fn[2])

    rEstimates <- resultA[,"estimate"]/resultB[,"estimate"]
    FNF <- fn/(fn+tp)
    rFNF <- FNF[1]/FNF[2] 
    FPF <- fp/(fp+tn)
    rFPF <-  FPF[1]/FPF[2] 
##    rEstimates <- c(rEstimates[1:2], rFNF, rFPF, rEstimates[3:4])
    rEstimates <- c(rEstimates[1:2], rEstimates[3:4])
    if(adjust)
        conf.level <- sqrt(conf.level)
    z <- qnorm((1 + conf.level) / 2)
##    EF <- exp(z*c(se.log.rTPF, se.log.rTNF, se.log.rFNF, se.log.rFPF, se.log.rLRp, se.log.rLRm))
    EF <- exp(z*c(se.log.rTPF, se.log.rTNF, se.log.rLRp, se.log.rLRm))

    res <- data.frame(matrix(NA, ncol=4, nrow=4))
    colnames(res) <- c("type", "lower", "estimate", "upper")
##     res[, 1] <- c("rTPF", "rTNF", "rFNF", "rFPF", "rLRplus", "rLRminus")
    res[, 1] <- c("rSens", "rSpec", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates/EF
    res[, 3] <- rEstimates
    res[, 4] <- rEstimates*EF
    res <- res[,c(1,3,2,4)]
    ## myorder <- c(1,4,2,3,5,6)
    ## res <- res[myorder,]
    return(res)

}
