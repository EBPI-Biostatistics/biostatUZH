#' Confidence intervals for the comparison of two diagnostic tests from
#' unpaired data
#' 
#' Compute confidence intervals for relative fractions and relative likelihood
#' ratios. Specifically, confidence intervals for the relative true positive,
#' true negative, false positive and false negative fraction as well as the
#' relative positive and negative likelihood ratio are provided.
#' 
#' 
#' @param tp Number of true positives of the two tests.
#' @param tn Number of true negatives.
#' @param fn Number of false negatives.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @param adjust Logical of length one indicating whether to compute adjusted
#' CIs or not. Default is FALSE.
#' @return A data.frame containing the estimated confidence intervals for the
#' measures of relative accuracy (first versus second test).
#' @author Leonhard Held
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford: Oxford University Press.
#' @examples
#' 
#' ## Calculate confidence intervals for data from a (hypothetical)
#' ## randomized unpaired study of chorionic villus sampling (CVS)
#' ## versus early amniocentesis (EA) for fetal abnormality from
#' ## Pepe (2003)
#' 
#' tp <- c(116, 111)
#' fp <- c(34, 111)
#' tn <- c(4844, 4765)
#' fn <- c(6, 13)
#' confIntIndependentDiagnostic(tp=tp, fp=fp, tn=tn, fn=fn)
#' 
confIntIndependentDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, adjust=FALSE)
{
    stopifnot(length(tp)==2, is.wholenumber(tp),
              length(fp)==2, is.wholenumber(fp),
              length(tn)==2, is.wholenumber(tn),
              length(fn)==2, is.wholenumber(fn),
              length(conf.level) ==1, 0 < conf.level, conf.level < 1,
              length(adjust)==1, is.logical(adjust))

    se.log.TPF <- sqrt(sum(1/tp) - sum(1/(tp+fn)))
    se.log.FPF <- sqrt(sum(1/fp) - sum(1/(fp+tn)))
    se.log.TNF <- sqrt(sum(1/tn) - sum(1/(tn+fp)))
    se.log.FNF <- sqrt(sum(1/fn) - sum(1/(fn+tp)))

    se.log.rTPF <- sqrt(sum(se.log.TPF^2))
    se.log.rTNF <- sqrt(sum(se.log.TNF^2))
    se.log.rFNF <- sqrt(sum(se.log.FNF^2))
    se.log.rFPF <- sqrt(sum(se.log.FPF^2))
    se.log.rLRp <- sqrt(sum(se.log.TPF^2 + se.log.FPF^2))
    se.log.rLRm <- sqrt(sum(se.log.TNF^2 + se.log.FNF^2))

    resultA <- confIntDiagnostic(tp[1], fp[1], tn[1], fn[1])
    resultB <- confIntDiagnostic(tp[2], fp[2], tn[2], fn[2])

    rEstimates <- resultA[,"estimate"] / resultB[,"estimate"]
    FNF <- fn / (fn + tp)
    rFNF <- FNF[1] / FNF[2] 
    FPF <- fp / (fp + tn)
    rFPF <-  FPF[1] / FPF[2] 
    rEstimates <- c(rEstimates[1:2], rEstimates[3:4])
    if(adjust)
        conf.level <- sqrt(conf.level)
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * c(se.log.rTPF, se.log.rTNF, se.log.rLRp, se.log.rLRm))

    res <- data.frame(matrix(NA, ncol=4, nrow=4))
    colnames(res) <- c("type", "lower", "estimate", "upper")
    res[, 1] <- c("rSens", "rSpec", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates/EF
    res[, 3] <- rEstimates
    res[, 4] <- rEstimates*EF
    res <- res[,c(1,3,2,4)]
    return(res)

}



#' Confidence intervals for operating characteristics of a diagnostic test
#' 
#' Compute confidence intervals for sensitivity, specificity, positive and
#' negative likelihood ratio and diagnostic odds ratio of a diagnostic test.
#' Optionally also positive and negative predictive value.
#' 
#' 
#' @param tp Number of true positives.
#' @param tn Number of true negatives.
#' @param fn Number of false negatives.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @param cohort Logical of length one indicating whether the data come from a
#' cohort study. Default is FALSE.
#' @param pr Prevalence.
#' @param digits Number of digits.
#' @return A dataframe containing the estimated confidence intervals for
#' sensitivity, specificity, positive and negative likelihood ratio. Optionally
#' also positive and negative predictive value.
#' @author Leonhard Held
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford: Oxford University Press.
#' @examples
#' 
#' ## Calculate confidence intervals for data from the Million Women Study
#' 
#' confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97)
#' confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, cohort=TRUE)
#' confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, pr=0.045)
#' confIntDiagnostic(tp=629, fp=3885, tn=117744, fn=97, digits=2)  
#' 
confIntDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, cohort=FALSE, pr=NA, digits=NA)
{
    stopifnot(is.wholenumber(tp), is.wholenumber(fp),
              is.wholenumber(tn), is.wholenumber(fn),  conf.level<1,
              conf.level>0)
    res <- data.frame(matrix(NA, nrow = 8, ncol = 4))

    colnames(res) <- c("type", "lower", "estimate", "upper")

    res[1, 2:4] <- wilson(x=tp, n=tp+fn, conf.level = conf.level)
    res[2, 2:4] <- wilson(x=tn, n=tn+fp, conf.level = conf.level)
    if((tp>0)&(fp>0)&(fn>0)&(tn>0)){
        LRplus <- confIntRiskRatio(x=c(tp,fp), n=c(tp+fn, fp+tn), conf.level = conf.level)
        LRminus <- confIntRiskRatio(x=c(fn,tn), n=c(tp+fn, tn+fp), conf.level = conf.level)
        DOR <- confIntOddsRatio(x=c(tp,fp), n=c(tp+fn, tn+fp), conf.level = conf.level)
        res[3, 2:4] <- LRplus
        res[4, 2:4] <- LRminus
        res[5, 2:4] <- DOR
    }
    
    if(!is.na(pr)){
        stopifnot(pr>0, pr<1)
        pr.odds <- pr/(1-pr)
        PPV.odds <- pr.odds*LRplus
        PPV <- PPV.odds/(1+PPV.odds)
        NPV.inv.odds <- pr.odds*LRminus
        NPV <- rev(1/(1+NPV.inv.odds))
        res[6, 2:4] <- PPV
        res[7, 2:4] <- NPV
        res[8, 3] <- pr
    }
    if(is.na(pr) & (cohort==TRUE)){
        res[6, 2:4] <- wilson(x=tp, n=tp+fp, conf.level = conf.level)
        res[7, 2:4] <- wilson(x=tn, n=tn+fn, conf.level = conf.level)
        res[8, 2:4] <- wilson(x=tp+fn, n=tp+tn+fp+fn, conf.level = conf.level)
    }

    res[, 1] <- c("Sensitivity", "Specificity", "LRplus", "LRminus", "DOR", "PPV", "NPV", "Prevalence")
    res <- res[,c(1,3,2,4)]
    
    if(!is.na(digits)){
        stopifnot(is.wholenumber(digits))
        res[, 2:4] <- as.numeric(sapply(res[, 2:4], 
                                        function(x) format(x, nsmall = digits, digits = digits)))
    }
    
    return(res)
    
}




#' Confidence intervals for the comparison of two diagnostic tests from paired
#' data
#' 
#' Compute confidence intervals for relative fractions. Specifically,
#' confidence intervals for the relative true positive, true negative, false
#' positive and false negative fraction are provided.
#' 
#' 
#' @param Diseased Frequency table with results from both tests in the diseased
#' population. First row and first column refers to negative test results.
#' @param nonDiseased Frequency table with results from both tests in the
#' non-diseased population. First row and first column refers to negative test
#' results.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @param adjust Logical of length one indicating whether to compute adjusted
#' CIs or not. Default is FALSE.
#' @return A data.frame containing the estimated confidence intervals for the
#' measures of relative accuracy (column test vs. row test).
#' @author Leonhard Held
#' @references Section 3.3 in Pepe, M.S. (2003) \emph{The statistical
#' evaluation of medical tests for classification and prediction}. Oxford
#' University Press.
#' @examples
#' 
#' ## Calculate confidence intervals for data from a study performing the
#' ## Exercise stress test (EST) and the determination of chest pain
#' ## history (CPH) for patients with suspected or probably coronary
#' ## heart disease (CHD).
#' 
#' Diseased <- matrix(c(25, 183, 29, 786), ncol = 2, nrow = 2, byrow = TRUE)
#' nonDiseased <- matrix(c(151, 176, 46, 69), ncol = 2, nrow = 2, byrow = TRUE)
#' colnames(Diseased) <- colnames(nonDiseased) <- c("CPH=0", "CPH=1")
#' rownames(Diseased) <- rownames(nonDiseased) <- c("EST=0", "EST=1")
#' 
#' confIntPairedDiagnostic(Diseased=Diseased, nonDiseased=nonDiseased)
#' 
confIntPairedDiagnostic <- function(Diseased, nonDiseased, conf.level = 0.95, adjust = FALSE)
{
    stopifnot(is.wholenumber(Diseased), is.wholenumber(nonDiseased),
              Diseased>0, nonDiseased>0,  conf.level<1,
              conf.level>0)
    stopifnot(nrow(Diseased)==2, ncol(Diseased)==2, nrow(nonDiseased)==2, ncol(nonDiseased)==2)

    rPF <- colSums(Diseased) / rowSums(Diseased)
    rNF <- colSums(nonDiseased) / rowSums(nonDiseased)
    rLR <- rPF/rNF
    
    se.log.rTPF <- sqrt((sum(Diseased)-sum(diag(Diseased)))/(rowSums(Diseased)[2]*colSums(Diseased)[2]))
    se.log.rFPF <- sqrt((sum(nonDiseased)-sum(diag(nonDiseased)))/(rowSums(nonDiseased)[2]*colSums(nonDiseased)[2]))
    se.log.rFNF <- sqrt((sum(Diseased)-sum(diag(Diseased)))/(rowSums(Diseased)[1]*colSums(Diseased)[1]))
    se.log.rTNF <- sqrt((sum(nonDiseased)-sum(diag(nonDiseased)))/(rowSums(nonDiseased)[1]*colSums(nonDiseased)[1]))
    se.log.rLRplus <- sqrt(se.log.rTPF^2 + se.log.rFPF^2)
    se.log.rLRminus <- sqrt(se.log.rFNF^2 + se.log.rTNF^2)
    
    rEstimates <- c(rPF[2], rNF[1], rev(rLR))
    if(adjust)
        conf.level <- sqrt(conf.level)
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z*c(se.log.rTPF, se.log.rTNF, se.log.rLRplus, se.log.rLRminus))

    res <- data.frame(matrix(NA, ncol=4, nrow=4))
    colnames(res) <- c("type", "lower", "estimate", "upper")
    res[, 1] <- c("rSens", "rSpec", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates/EF
    res[, 3] <- rEstimates
    res[, 4] <- rEstimates*EF
    res <- res[,c(1,3,2,4)]

    return(res)
}
