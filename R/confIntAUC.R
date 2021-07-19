#' Confidence interval for AUC
#' 
#' Compute Wald confidence intervals for the area under the curve on the
#' original and logit scale.
#' 
#' 
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A list with the entries: \describe{ \item{p}{Estimated proportion.}
#' \item{CIs}{data.frame containing the estimated confidence intervals.} }
#' @author Leonhard Held
#' @references Altman, D.G., Machin, D., Bryant, T.N. and Gardner, M.J. (2001).
#' \emph{Statistics with confidence}. 2nd Edition, 2000. BMJ Books. Chapter 10.
#' @examples
#' 
#' set.seed(12345)
#' cases <- rnorm(100, mean=2)
#' controls <- rnorm(50)
#' confIntAUC(cases, controls)
#' 
confIntAUC <- function(cases, controls, conf.level = 0.95){
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






#' Confidence interval for the difference in AUC based on two independent
#' samples
#' 
#' Computes confidence interval for the difference in the area under the curve
#' based on two independent samples.
#' 
#' For type="Wald", standard Wald confidence intervals are calculated for AUC
#' of both tests and their difference. For type="logit", the substitution
#' method is used based on the logit transformation for the AUC of both tests.
#' The confidence interval for the difference in AUC is then calculated using
#' Newcombe's method.
#' 
#' @param casesA Values of the continuous variable from Test A for the cases.
#' @param controlsA Values of the continuous variable from Test A for the
#' controls.
#' @param casesB Values of the continuous variable from Test B for the cases.
#' @param controlsB Values of the continuous variable from Test B for the
#' controls.
#' @param type "Wald" (default) or "Logit".
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A matrix with estimate and confidence limits for AUC of the two
#' tests and their difference.
#' @author Leonhard Held
#' @references Newcombe, R.G. (1998). \emph{Interval estimation for the
#' difference between independent proportions: Comparison of eleven methods.}
#' Stat. Med., *17*, 873-890.
#' 
#' Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for
#' classification and prediction}. Oxford University Press.
#' @examples
#' 
#' set.seed(12345)
#' casesA <- rnorm(200, mean=2.5)
#' controlsA <- rnorm(100)
#' casesB <- rnorm(100, mean=1.5)
#' controlsB <- rnorm(200)
#' 
#' confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type="Wald")
#' confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type="Logit")
#' 
confIntIndependentAUCDiff <- function(casesA, controlsA, casesB, controlsB,
                                      type="Wald", conf.level = 0.95)
{
    ncontrolsA <- nrow(controlsA)
    ncasesA <- nrow(casesA)
    ncontrolsB <- nrow(controlsB)
    ncasesB <- nrow(casesB)
    auc <- numeric()

    resA <- confIntAUC(casesA, controlsA, conf.level = conf.level)
    resB <- confIntAUC(casesB, controlsB, conf.level = conf.level)
    factor <- qnorm((1 + conf.level) / 2)

    type <- match.arg(type, c("Wald", "Logit"))

    if(type=="Wald"){
        ## take intervals on original scale
        aucA <- resA[1,3]
        aucB <- resB[1,3]
        D <- aucA - aucB
        lowerA <- resA[1,2]
        lowerB <- resB[1,2]
        upperA <- resA[1,4]
        upperB <- resB[1,4]
        ## compute and combine standard errors
        seA <- (upperA-lowerA)/(2*factor)
        seB <- (upperB-lowerB)/(2*factor)
        se.D <- sqrt(seA^2+seB^2)
        D.lower <- D - factor * se.D
        D.upper <- D + factor * se.D

        res <- data.frame(matrix(NA, ncol = 4))
        colnames(res) <- c("outcome", "lower", "estimate", "upper")
        res[1, 2:4] <- resA[1,2:4] # Wald interval on original scale
        res[2, 2:4] <- resB[1,2:4] # Wald interval on original scale
        res[3, 2:4] <- c(D.lower, D, D.upper)
        res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    }

    if(type=="Logit"){
        ## take intervals on logit scale
        aucA <- resA[2,3]
        aucB <- resB[2,3]
        D <- aucA - aucB
        lowerA <- resA[2,2]
        lowerB <- resB[2,2]
        upperA <- resA[2,4]
        upperB <- resB[2,4]

        ## Apply Newcombe trick
        D.lower <- D - sqrt((aucA - lowerA) ^ 2 + (aucB - upperB) ^ 2)
        D.upper <- D + sqrt((aucA - upperA) ^ 2 + (aucB - lowerB) ^ 2)

        res <- data.frame(matrix(NA, ncol = 4))
        colnames(res) <- c("outcome", "lower", "estimate", "upper")
        res[1, 2:4] <- resA[2,2:4] # Wald interval on logit scale
        res[2, 2:4] <- resB[2,2:4] # Wald interval on logit scale
        res[3, 2:4] <- c(D.lower, D, D.upper)
        res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    }

    return(res)
}



#' Confidence interval for the difference in AUC based on two paired samples
#' 
#' Computes confidence interval for the difference in the area under the curve
#' based on two paired samples.
#' 
#' 
#' @param cases Matrix with values of the continuous variable for the cases.
#' First column gives values of test A, second gives values of test B.
#' @param controls Matrix with values of the continuous variable for the
#' controls. First column gives values of test A, second gives values of test
#' B.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A matrix with estimate and confidence limits for AUC of the two
#' tests and their difference.
#' @author Leonhard Held
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford University Press.
#' @examples
#' 
#' data(wiedat2b)
#' ind <- wiedat2b[,3]
#' cases <- wiedat2b[ind==1, 1:2]
#' controls <- wiedat2b[ind==0, 1:2]
#' confIntPairedAUCDiff(cases, controls)
#' 
confIntPairedAUCDiff <- function(cases, controls, conf.level = 0.95){
    stopifnot(is.matrix(cases) || is.data.frame(cases),
              is.matrix(controls) || is.data.frame(controls))

    # estimate AUC as normalized test statistic of Wilcoxon test
    ncontrols <- nrow(controls)
    ncases <- nrow(cases)
    auc <- numeric()
    for(k in 1:2)
        auc[k] <- as.numeric(wilcox.test(cases[,k], controls[,k], exact = FALSE)$statistic / (ncases * ncontrols))
    auc.diff <- auc[1] - auc[2]
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
