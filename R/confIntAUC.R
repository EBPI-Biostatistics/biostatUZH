#' Confidence interval for AUC
#' 
#' Compute Wald confidence intervals for the area under the curve on the
#' original and logit scale.
#' 
#' 
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A data.frame with the columns:
#' \describe{
#' \item{type:}{Type of confidence interval.}
#' \item{lower:}{Lower bound of CI.}
#' \item{AUC:}{Area under the curve.}
#' \item{upper}{Upper bound of CI.}}
#' @author Leonhard Held
#' @seealso \code{\link{confIntIndependentAUCDiff}},
#' \code{\link{confIntPairedAUCDiff}}, \code{\link{standardErrorAUCDiff}}
#' @references Altman, D.G., Machin, D., Bryant, T.N. and Gardner, M.J. (2001).
#' \emph{Statistics with confidence}. 2nd Edition, 2000. BMJ Books. Chapter 10.
#' @examples
#' 
#' set.seed(12345)
#' cases <- rnorm(100, mean=2)
#' controls <- rnorm(50)
#' confIntAUC(cases, controls)
#' 
#' @export
confIntAUC <- function(cases, controls, conf.level = 0.95){

    stopifnot(is.numeric(cases), length(cases) >= 1, is.finite(cases),
              is.numeric(controls), length(controls) >= 1, is.finite(controls),
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1)

    ## estimate AUC as normalized test statistic of Wilcoxon test
    ncontrols <- length(controls)
    ncases <- length(cases)
    auc <- as.numeric(wilcox.test(cases, controls, exact = FALSE)$statistic /
               (ncases * ncontrols))
    aucSE <- standardErrorAUC(cases, controls)

    ## compute confidence intervals on original scale
    z <- qnorm((1 + conf.level) / 2)
    lower <- auc - z * aucSE
    upper <- auc + z * aucSE

    ## on logit scale
    logitAuc <- log(auc / (1 - auc))
    logitAucSE <- aucSE / (auc * (1 - auc))
    logitLowCI <- logitAuc - z * logitAucSE
    logitUpCI <- logitAuc + z * logitAucSE

    ## backtransformation
    lowerLogit <- 1 / (1 + exp(-logitLowCI))
    upperLogit <- 1 / (1 + exp(-logitUpCI))
    
    res <- data.frame(matrix(NA, ncol = 4))
    colnames(res) <- c("type", "lower", "AUC", "upper")
    res[1, 2:4] <- c(lower, auc, upper)
    res[2, 2:4] <- c(lowerLogit, auc, upperLogit)
    res[, 1] <- c("Wald", "logit Wald")
    res
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
#' @param type "Wald" (default) or "logit".
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @return A data.frame with estimate and confidence limits for AUC of the two
#' tests and their difference.
#' @author Leonhard Held
#' @seealso \code{\link{confIntAUC}}, \code{\link{confIntPairedAUCDiff}},
#' \code{\link{standardErrorAUCDiff}}
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
#' confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type = "Wald")
#' confIntIndependentAUCDiff(casesA, controlsA, casesB, controlsB, type = "logit")
#' 
#' @export
confIntIndependentAUCDiff <- function(casesA, controlsA, casesB, controlsB,
                                      type = c("Wald", "logit"), conf.level = 0.95)
{
    stopifnot(is.numeric(casesA), length(casesA) >= 1, is.finite(casesA),
              is.numeric(controlsA), length(controlsA) >= 1, is.finite(controlsA),
              is.numeric(casesB), length(casesB) >= 1, is.finite(casesB),
              is.numeric(controlsB), length(controlsB) >= 1, is.finite(controlsB),
              !is.null(type))
    type <- match.arg(tolower(type), choices = c("wald", "logit"))
    stopifnot(is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1)
    
    resA <- confIntAUC(casesA, controlsA, conf.level = conf.level)
    resB <- confIntAUC(casesB, controlsB, conf.level = conf.level)
    factor <- qnorm((1 + conf.level) / 2)

    if(type == "Wald"){
        ## take intervals on original scale
        aucA <- resA[1, 3]
        aucB <- resB[1, 3]
        D <- aucA - aucB
        lowerA <- resA[1, 2]
        lowerB <- resB[1, 2]
        upperA <- resA[1, 4]
        upperB <- resB[1, 4]
        ## compute and combine standard errors
        seA <- (upperA - lowerA) / (2 * factor)
        seB <- (upperB - lowerB) / (2 * factor)
        seD <- sqrt(seA^2 + seB^2)
        lowerD <- D - factor * seD
        upperD <- D + factor * seD

        res <- data.frame(matrix(data = NA, ncol = 4))
        colnames(res) <- c("outcome", "lower", "estimate", "upper")
        res[1, 2:4] <- resA[1, 2:4] # Wald interval on original scale
        res[2, 2:4] <- resB[1, 2:4] # Wald interval on original scale
        res[3, 2:4] <- c(lowerD, D, upperD)
        res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
        return(res)
    }
    ## else type == "Logit"
    ## take intervals on logit scale
    aucA <- resA[2, 3]
    aucB <- resB[2, 3]
    D <- aucA - aucB
    lowerA <- resA[2, 2]
    lowerB <- resB[2, 2]
    upperA <- resA[2, 4]
    upperB <- resB[2, 4]
    
    ## Apply Newcombe trick
    lowerD <- D - sqrt((aucA - lowerA)^2 + (aucB - upperB)^2)
    upperD <- D + sqrt((aucA - upperA)^2 + (aucB - lowerB)^2)
    
    res <- data.frame(matrix(data = NA, ncol = 4))
    colnames(res) <- c("outcome", "lower", "estimate", "upper")
    res[1, 2:4] <- resA[2, 2:4] # Wald interval on logit scale
    res[2, 2:4] <- resB[2, 2:4] # Wald interval on logit scale
    res[3, 2:4] <- c(lowerD, D, upperD)
    res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    res
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
#' @return A data.frame with the estimates and confidence limits for the AUC of the two
#' tests and their difference.
#' @author Leonhard Held
#' @seealso \code{\link{confIntAUC}}, \code{\link{confIntIndependentAUCDiff}},
#' \code{\link{standardErrorAUCDiff}}
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford University Press.
#' @examples
#' 
#' data(wiedat2b)
#' ind <- wiedat2b[,3]
#' cases <- wiedat2b[ind == 1, 1:2]
#' controls <- wiedat2b[ind == 0, 1:2]
#' confIntPairedAUCDiff(cases = cases, controls = controls)
#' 
#' @export
confIntPairedAUCDiff <- function(cases, controls, conf.level = 0.95){
    cases <- as.matrix(cases[,1:2])
    stopifnot(is.numeric(cases), length(cases) > 1,
              ncol(cases) == 2, is.finite(cases))
    controls <- as.matrix(controls[,1:2])
    stopifnot(is.numeric(controls), length(controls) > 1,
              ncol(controls) >= 2, is.finite(controls),
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1)
    
    # estimate AUC as normalized test statistic of Wilcoxon test
    ncontrols <- nrow(controls)
    ncases <- nrow(cases)
    auc <- numeric(2)
    for(k in 1:2)
        auc[k] <- as.numeric(wilcox.test(cases[,k], controls[,k], exact = FALSE)$statistic / (ncases * ncontrols))
    aucDiff <- auc[1] - auc[2]
    aucDiffSE <- standardErrorAUCDiff(cases, controls)

    # compute confidence intervals on original scale
    z <- qnorm((1 + conf.level) / 2)
    lower <- aucDiff - z * aucDiffSE
    upper <- aucDiff + z * aucDiffSE

    res <- data.frame(matrix(NA, ncol = 4))
    colnames(res) <- c("outcome", "lower", "estimate", "upper")
    res[1, 2:4] <- confIntAUC(cases[,1], controls[,1])[2, 2:4] ## avoids overshoot
    res[2, 2:4] <- confIntAUC(cases[,2], controls[,2])[2, 2:4] ## avoids overshoot
    res[3, 2:4] <- c(lower, aucDiff, upper)
    res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    res
}


#' Standard Error of a AUC difference
#' 
#' Computes the standard error of the difference in area under the curve for
#' paired samples.
#' 
#' 
#' @param cases Matrix with values of the continuous variable for the cases.
#' First column gives values of test A, second gives values of test B.
#' @param controls Matrix with values of the continuous variable for the
#' controls. First column gives values of test A, second gives values of test B.
#' @return The standard error.
#' @author Leonhard Held
#' @seealso \code{\link{confIntAUC}}, \code{\link{confIntIndependentAUCDiff}},
#' \code{\link{confIntPairedAUCDiff}}
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford University Press.
#' @keywords univar htest
#' @export
standardErrorAUCDiff <- function(cases, controls){
    cases <- as.matrix(cases[,1:2])
    stopifnot(is.numeric(cases), length(cases) > 1,
              ncol(cases) == 2, is.finite(cases))
    controls <- as.matrix(controls[,1:2])
    stopifnot(is.numeric(controls), length(controls) > 1,
              ncol(controls) >= 2, is.finite(controls))

    ncases <- nrow(cases)
    ncontrols <- nrow(controls)
    # non-disease placement values of cases
    C <- matrix(data = NA, nrow = ncases, ncol = 2)
    # disease placement values of controls
    R <- matrix(data = NA, nrow = ncontrols, ncol = 2)

    for(k in 1:2){
        for(i in 1:ncases)
            C[i,k] <- mean(as.numeric(controls[,k] < cases[i,k]) + 0.5 * as.numeric(controls[,k] == cases[i,k]))
        for(j in 1:ncontrols)
            R[j,k] <- mean(as.numeric(cases[,k] > controls[j,k]) + 0.5 * as.numeric(cases[,k] == controls[j,k]))
    }
    sqrt((var(R[,1] - R[,2]) / ncontrols + var(C[,1] - C[,2]) / ncases))
}
