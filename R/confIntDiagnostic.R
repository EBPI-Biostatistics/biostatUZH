#' Confidence intervals for the comparison of two diagnostic tests from
#' unpaired data
#' 
#' Compute confidence intervals for relative fractions and relative likelihood
#' ratios. Specifically, confidence intervals for the relative true positive,
#' true negative, false positive and false negative fraction as well as the
#' relative positive and negative likelihood ratio are provided.
#' 
#' 
#' @param tp Vector of length 2. Number of true positives of the two tests.
#' @param fp Vector of length 2. Number of false negatives of the two tests.
#' @param tn Vector of length 2. Number of true negatives of the two tests.
#' @param fn Vector of length 2. Number of false negatives of the two tests.
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
#' confIntIndependentDiagnostic(tp = tp, fp = fp, tn = tn, fn = fn)
#' 
#' @export
confIntIndependentDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, adjust=FALSE)
{
    stopifnot(is.numeric(tp), length(tp) == 2, is.finite(tp), is.wholenumber(tp),
              is.numeric(fp), length(fp) == 2, is.finite(fp), is.wholenumber(fp),
              is.numeric(tn), length(tn) == 2, is.finite(tn), is.wholenumber(tn),
              is.numeric(fn), length(fn) == 2, is.finite(fn), is.wholenumber(fn),
              is.numeric(conf.level), length(conf.level) == 1, is.finite(conf.level), 
              0 < conf.level, conf.level < 1,
              is.logical(adjust), length(adjust) == 1, is.finite(adjust))

    se.log.TPF <- sqrt(sum(1 / tp) - sum(1 / (tp + fn)))
    se.log.FPF <- sqrt(sum(1 / fp) - sum(1 / (fp + tn)))
    se.log.TNF <- sqrt(sum(1 / tn) - sum(1 / (tn + fp)))
    se.log.FNF <- sqrt(sum(1 / fn) - sum(1 / (fn + tp)))

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
    colnames(res) <- c("type", "estimate", "lower", "upper")
    res[, 1] <- c("rSens", "rSpec", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates
    res[, 3] <- rEstimates / EF
    res[, 4] <- rEstimates * EF
    res
}



#' Confidence intervals for operating characteristics of a diagnostic test
#' 
#' Compute confidence intervals for sensitivity, specificity, positive and
#' negative likelihood ratio and diagnostic odds ratio of a diagnostic test.
#' Optionally also positive and negative predictive value.
#' 
#' 
#' @param tp Number of true positives.
#' @param fp Number of false negatives.
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
#' confIntDiagnostic(tp = 629, fp = 3885, tn = 117744, fn = 97)
#' confIntDiagnostic(tp = 629, fp = 3885, tn = 117744, fn = 97, cohort = TRUE)
#' confIntDiagnostic(tp = 629, fp = 3885, tn = 117744, fn = 97, pr = 0.045)
#' confIntDiagnostic(tp = 629, fp = 3885, tn = 117744, fn = 97, digits = 2)  
#' 
#' @export
confIntDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, cohort = FALSE, pr = NA,
                              digits = NA)
{
    stopifnot(is.numeric(tp), length(tp) == 1, is.finite(tp), is.wholenumber(tp),
              is.numeric(fp), length(fp) == 1, is.finite(fp), is.wholenumber(fp),
              is.numeric(tn), length(tn) == 1, is.finite(tn), is.wholenumber(tn),
              is.numeric(fn), length(fn) == 1, is.finite(fn), is.wholenumber(fn),
              is.numeric(conf.level), length(conf.level) == 1, is.finite(conf.level), 
              0 < conf.level, conf.level < 1,
              is.logical(cohort), length(cohort)==1, is.finite(cohort),
              length(pr) == 1, is.na(pr) || (0 < pr && pr < 1),
              length(digits) == 1,
              is.na(digits) || (is.numeric(digits) && is.finite(digits) &&
                                is.wholenumber(digits) && digits >= 0)
              )

    res <- data.frame(matrix(data = NA, nrow = 8, ncol = 4))

    colnames(res) <- c("type", "lower", "estimate", "upper")

    res[1, 2:4] <- wilson(x = tp, n = tp + fn, conf.level = conf.level)
    res[2, 2:4] <- wilson(x = tn, n= tn + fp, conf.level = conf.level)
    if(tp > 0 && fp > 0 && fn > 0 && tn > 0){
        LRplus <- confIntRiskRatio(x = c(tp, fp), n = c(tp + fn, fp + tn),
                                   conf.level = conf.level)
        LRminus <- confIntRiskRatio(x = c(fn, tn), n = c(tp + fn, tn + fp),
                                    conf.level = conf.level)
        DOR <- confIntOddsRatio(x = c(tp, fp), n = c(tp + fn, tn + fp),
                                conf.level = conf.level)
        res[3, 2:4] <- LRplus 
        res[4, 2:4] <- LRminus
        res[5, 2:4] <- DOR
    }
    
    if(!is.na(pr)){
        pr.odds <- pr / (1 - pr)
        PPV.odds <- pr.odds * LRplus
        PPV <- PPV.odds/(1 + PPV.odds)
        NPV.inv.odds <- pr.odds * LRminus
        NPV <- rev(1 / (1 + NPV.inv.odds))
        res[6, 2:4] <- PPV
        res[7, 2:4] <- NPV
        res[8, 3] <- pr
    }
    if(is.na(pr) && cohort){
        res[6, 2:4] <- wilson(x = tp, n = tp + fp, conf.level = conf.level)
        res[7, 2:4] <- wilson(x = tn, n = tn + fn, conf.level = conf.level)
        res[8, 2:4] <- wilson(x = tp + fn, n = tp + tn + fp + fn, conf.level = conf.level)
    }

    res[, 1] <- c("Sensitivity", "Specificity", "LRplus", "LRminus", "DOR", "PPV",
                  "NPV", "Prevalence")
    res <- res[, c(1,3,2,4)]
    
    if(!is.na(digits)){
        formatedStr <- format(unlist(res[, 2:4]), nsmall = digits, digits = digits)
        res[, 2:4] <- as.numeric(ifelse(is.na(unlist(res[, 2:4])), NA, formatedStr))
    }
    res
}




#' Confidence intervals for the comparison of two diagnostic tests from paired
#' data
#' 
#' Compute confidence intervals for relative fractions. Specifically,
#' confidence intervals for the relative true positive, true negative, false
#' positive and false negative fraction are provided.
#' 
#' 
#' @param diseased 2 by 2 frequency table with results from both tests in the diseased
#' population. First row and first column refers to negative test results.
#' @param nonDiseased 2 by 2 frequency table with results from both tests in the
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
#' diseased <- matrix(c(25, 183, 29, 786), ncol = 2, nrow = 2, byrow = TRUE)
#' nonDiseased <- matrix(c(151, 176, 46, 69), ncol = 2, nrow = 2, byrow = TRUE)
#' colnames(diseased) <- colnames(nonDiseased) <- c("CPH=0", "CPH=1")
#' rownames(diseased) <- rownames(nonDiseased) <- c("EST=0", "EST=1")
#' 
#' confIntPairedDiagnostic(diseased = diseased, nonDiseased = nonDiseased)
#' 
#' @export
confIntPairedDiagnostic <- function(diseased, nonDiseased, conf.level = 0.95, adjust = FALSE)
{
    stopifnot(is.numeric(diseased), is.matrix(diseased), dim(diseased) == c(2, 2),
              is.finite(diseased), is.wholenumber(diseased), diseased > 0,
              is.numeric(nonDiseased), is.matrix(nonDiseased), dim(nonDiseased) == c(2, 2),
              is.finite(nonDiseased), is.wholenumber(nonDiseased), nonDiseased > 0,
              is.numeric(conf.level), length(conf.level) == 1, is.finite(conf.level), 
              0 < conf.level, conf.level < 1,
              is.logical(adjust), length(adjust) == 1, is.finite(adjust))

    rPF <- colSums(diseased) / rowSums(diseased)
    rNF <- colSums(nonDiseased) / rowSums(nonDiseased)
    rLR <- rPF / rNF
    
    seLogrTPF <- sqrt((sum(diseased) - sum(diag(diseased))) /
                        (rowSums(diseased)[2] * colSums(diseased)[2]))
    seLogrFPF <- sqrt((sum(nonDiseased) - sum(diag(nonDiseased))) /
                        (rowSums(nonDiseased)[2] * colSums(nonDiseased)[2]))
    seLogrFNF <- sqrt((sum(diseased) - sum(diag(diseased))) /
                        (rowSums(diseased)[1] * colSums(diseased)[1]))
    seLogrTNF <- sqrt((sum(nonDiseased) - sum(diag(nonDiseased))) /
                        (rowSums(nonDiseased)[1] * colSums(nonDiseased)[1]))
    seLogrLRplus <- sqrt(seLogrTPF^2 + seLogrFPF^2)
    seLogrLRminus <- sqrt(seLogrFNF^2 + seLogrTNF^2)
    
    rEstimates <- c(rPF[2], rNF[1], rev(rLR))
    if(adjust)
        conf.level <- sqrt(conf.level)
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z*c(seLogrTPF, seLogrTNF, seLogrLRplus, seLogrLRminus))

    res <- data.frame(matrix(data = NA, ncol = 4, nrow = 4))
    colnames(res) <- c("type", "estimate", "lower", "upper")
    res[, 1] <- c("rSens", "rSpec", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates
    res[, 3] <- rEstimates / EF
    res[, 4] <- rEstimates * EF
    res
}


#' Confidence interval for a risk difference
#' 
#' Compute a confidence interval for a risk difference based on
#' Wald and Wilson confidence intervals for the individual risks.
#' 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A list with the entries: \itemize{ \item rd: Estimated risk difference.
#' \item CIs: Data.frame containing confidence intervals for the risk difference. }
#' @author Leonhard Held
#' @seealso \code{\link{wilson}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntOddsRatio}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 6
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntRiskDiff(x = x, n = n)$CIs
#' 
#' @export
confIntRiskDiff <- function(x, n, conf.level = 0.95){
    stopifnot(is.numeric(x), length(x) == 2,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 2,
              is.finite(n), is.wholenumber(n),
              0 < x, x < n,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    ci1 <-  wilson(x[1], n[1], conf.level = conf.level)
    ci2 <-  wilson(x[2], n[2], conf.level = conf.level)

    diff <- matrix(ci1[2] - ci2[2])
    se.diff <- sqrt((ci1[2] * (1 - ci1[2])) / n[1] + (ci2[2] * (1 - ci2[2])) / n[2])  
    z <- qnorm((1 + conf.level) / 2)
    wald.lower <- diff - z * se.diff
    wald.upper <- diff + z * se.diff
    score.lower <- diff - sqrt((ci1[2] - ci1[1])^2 + (ci2[3] - ci2[2])^2)
    score.upper <- diff + sqrt((ci2[2] - ci2[1])^2 + (ci1[3] - ci1[2])^2)

    result <- matrix(ncol=2, nrow=2)
    result[,1] <- c(wald.lower, score.lower)
    result[,2] <- c(wald.upper, score.upper)

    out <- data.frame(type = c("Wald", "Wilson"), result)
    names(out) <- c("type", "lower", "upper")

    list("rd" = diff, "CIs" = out)
}




#' Confidence interval for a risk ratio
#' 
#' Provides a confidence interval for a risk ratio. The method is based on the
#' Wald interval for the log risk ratio.
#' 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A vector containing the risk ratio and the limits of the confidence
#' interval.
#' @author Leonhard Held
#' @seealso \code{\link{confIntDiagnostic}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntOddsRatio}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 7
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntRiskRatio(x = x, n = n)
#' 
#' @export
confIntRiskRatio <- function(x, n, conf.level = 0.95){
    stopifnot(is.numeric(x), length(x) == 2,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 2,
              is.finite(n), is.wholenumber(n),
              0 < x, x < n,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    Risk <- x / n
    RiskRatio <- Risk[1] / Risk[2]
    se.log.RiskRatio <- sqrt(sum(1 / x) - sum(1 / n))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * se.log.RiskRatio)
    wald.lower <- RiskRatio / EF
    wald.upper <- RiskRatio * EF

    c("lower" = wald.lower, "Risk Ratio" = RiskRatio, "upper" = wald.upper)
}




#' Confidence interval for an odds ratio
#' 
#' Confidence interval for an odds ratio. The method is based on the
#' Wald interval for the log odds ratio. 
#' 
#' @param x Vector of length 2, number of successes in each group.
#' @param n Vector of length 2, total number of trials in each group.
#' @param conf.level Confidence level for confidence interval. Default value is
#' 0.95.
#' @return A vector containing the odds ratio and the limits of the confidence
#' interval.
#' @author Leonhard Held
#' @seealso \code{\link{confIntDiagnostic}}, \code{\link{confIntRiskRatio}},
#' \code{\link{confIntRiskDiff}}.
#' @references DG Altman, D Machin, TN Bryant, MJ Gardner. \emph{Statistics
#' with confidence}, 2nd Edition, 2000, Chapter 7
#' @keywords htest
#' @examples
#' 
#' x <- c(30, 50)
#' n <- c(100, 120)
#' confIntOddsRatio(x = x, n = n)
#' 
#' @export
confIntOddsRatio <- function(x, n, conf.level = 0.95){
    stopifnot(is.numeric(x), length(x) == 2,
              is.finite(x), is.wholenumber(x),
              is.numeric(n), length(n) == 2,
              is.finite(n), is.wholenumber(n),
              0 < x, x < n,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    y <- n - x
    Odds <- x / y
    OddsRatio <- Odds[1] / Odds[2]
    se.log.OddsRatio <- sqrt(sum(1 / x) + sum(1 / y))
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z * se.log.OddsRatio)
    wald.lower <- OddsRatio / EF
    wald.upper <- OddsRatio * EF

    c("lower" = wald.lower, "Odds Ratio" = OddsRatio, "upper" = wald.upper)
}
