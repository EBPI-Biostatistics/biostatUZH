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



#' ROC curve and (bootstrap) confidence interval for AUC in the binormal model
#' 
#' This function computes for values of a continuous variable of a group of
#' cases and a group of controls the ROC curve in the binormal model.
#' Additionally, the AUC including a confidence interval is provided, with
#' bootstrap or with Wald, assuming the variances are equal (the assumption of
#' unequal variances has yet to be implemented).
#' 
#' @aliases confIntAUCbinorm AUCbinorm
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param conf.level Confidence level for confidence interval.
#' @param method Either "boot" for bootstrap CI or "wald" for a Wald CI.
#' @param replicates Number of boostrap replicates. Only used if
#' \code{method = 'boot'}.
#' @param grid Number of intervals to split \eqn{[0, 1]}.
#' Used to compute the ROC curve if \code{method = 'boot'}.
#' @param var.equal Logical with default TRUE.
#' Are the variances assumed to be equal or not?
#' Only used if \code{method = 'wald'}.
#' Currently, only \code{var.equal = TRUE} is supported.
#' @return The results for \code{method = 'boot'} are
#' \item{a}{The values of a.}
#' \item{b}{The values of b.}
#' \item{x.val}{The values on the \eqn{x}-axis of a ROC plot.}
#' \item{y.val}{The values on the \eqn{y}-axis of a ROC plot, i.e. the binormal ROC curve.}
#' \item{auc}{Area under the ROC curve.}
#' \item{lowBootCI}{Lower limit of bootstrap confidence interval.}
#' \item{upBootCI}{Upper limit of bootstrap confidence interval.}
#' 
#' The results for \code{method = 'wald'} are
#' \item{a}{The values of a.}
#' \item{b}{The values of b.}
#' \item{auc}{Area under the ROC curve.}
#' \item{lowCI}{Lower limit of confidence interval.}
#' \item{upCI}{Upper limit of confidence interval.}
#' @note The Wald approach is documented in a vignette: see
#' \code{vignette("aucbinormal")}
#' @author Kaspar Rufibach (\code{method = 'boot'}) and Leonhard Held
#' (\code{method = 'wald'})
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford: Oxford University Press.
#' @seealso \code{\link{summaryROC}}
#' @keywords htest
#' @examples
#'
#' ## simulate data
#' ## --------------
#' set.seed(1977)
#' controls <- rnorm(n = 50)
#' cases <-  rnorm(n = 40, mean = 0.5, sd = 1.5)
#' 
#' ## summary of ROC curve
#' ## --------------
#' res <- summaryROC(cases, controls, conf.level = 0.95)
#' 
#' 
#' ## alternative bootstrap CI for AUC
#' ## --------------
#' resBinormBoot <- confIntAUCbinorm(cases = cases, controls = controls, 
#'                                   conf.level = 0.95, replicates = 1000,
#'                                   grid = 100)
#' 
#' ## alternative bootstrap CI for AUC 
#' ## --------------
#' resBinormWald <- confIntAUCbinorm(cases, controls, conf.level = 0.95,
#'                                   method = "wald")
#' 
#' ## display results
#' ## --------------
#' str(res)
#' str(resBinormBoot)
#' str(resBinormWald)
#' 
#' ## plot ROC curve
#' ## --------------
#' plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
#'      xlab = "1 - specificity", ylab = "sensitivity", pty = 's')
#' segments(0, 0, 1, 1, lty = 2)
#' lines(res$x.val, res$y.val, type = 'l', col = 2, lwd = 2, lty = 2)
#' lines(resBinormBoot$x.val, resBinormBoot$y.val, type = 'l',
#'       col = 4, lwd = 2, lty = 2)
#' 
#' @importFrom boot boot boot.ci
#' @export
confIntAUCbinorm <- function(cases, controls, conf.level = 0.95, method = c("boot", "wald"),
                             replicates = 1000, grid = 100, var.equal = TRUE){ 

    stopifnot(is.numeric(cases), length(cases) >= 1, is.finite(cases),
              is.numeric(controls), length(controls) >= 1, is.finite(controls),
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1,
              !is.null(method))
    method <- match.arg(method)
    stopifnot(is.numeric(replicates), length(replicates) == 1,
              is.wholenumber(replicates), replicates >= 1,
              is.numeric(grid), length(grid) == 1,
              is.wholenumber(grid), grid >= 1,
              is.logical(var.equal), is.finite(var.equal),
              length(var.equal) == 1)
    
    alpha <- 1 - conf.level

    ## Bootstrap CI
    ## -------------------
    if (method == "boot") {
        ## compute bootstrap CI for binormal AUC
        ## data in two columns: y (measurements), status (0 = control, 1 = case)
        
        ts <- seq(from = 0, to = 1, length.out = grid)		
        a <- (mean(cases) - mean(controls)) / sd(cases)
        b <- sd(controls) / sd(cases)	
        biroc <- pnorm(a + b * qnorm(ts))

        dat <- cbind(c(controls, cases), c(rep(0, length(controls)), rep(1, length(cases))))

        AUCbinorm <- function(data, x) {
            dat <- data[x, ]
            cases <- dat[dat[, 2] == 1, 1]
            controls <- dat[dat[, 2] == 0, 1]
            a <- (mean(cases) - mean(controls)) / sd(cases)
            b <- sd(controls) / sd(cases)
            auc <- pnorm(a * (b ^ 2 + 1) ^ (-1 / 2))
            auc
        }

        ## compute bootstrap ci
        res <- boot(data = dat, statistic = AUCbinorm, R = replicates)
        aucbin <- boot.ci(boot.out = res, type = "bca")

        return(list("a" = a, "b" = b, "x.val" = ts, "y.val" = biroc,
                    "auc" = res$t0, "lowCI" = aucbin$bca[4],
                    "upCI" = aucbin$bca[5], method = method))  
    }

    ## Wald CI
    ## -------------------
    if(method == "wald") {
        if(!var.equal)
            stop("var.equal = FALSE not implemented for method = 'wald'")
        n0 <- length(controls)
        n1 <- length(cases)
        mu0 <- mean(controls)
        mu1 <- mean(cases)
        
        s0 <- s1 <- sd(c(cases, controls))
        auc.se <- sqrt(1 / n0 + 1 / n1)
        
        a <- (mu1 - mu0) / s1
        b <- s0 / s1
        auc <- pnorm(a / (sqrt(1 + b^2)))
        lowCI <- pnorm((a - qnorm(1 - alpha / 2) * auc.se) / (sqrt(1 + b^2)))
        upCI <- pnorm((a + qnorm(1-alpha/2) * auc.se) / (sqrt(1 + b^2)))
        
        return(list("a" = a, "b" = b, "auc" = auc, "lowCI" = lowCI, "upCI" = upCI,
                    method = method))
    }
}


#' ROC curve and an asymptotic confidence interval for AUC
#' 
#' This function computes the ROC curve for values of a continuous variable
#' of a group of cases and a group of controls the. Additionally, the AUC
#' with an asymptotic confidence interval is provided.
#'
#' 
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param conf.level Confidence level of confidence interval.
#' @return
#' \item{x.val}{1-specificity of the test, so the values on the \eqn{x}-axis
#' of a ROC plot.}
#' \item{y.val}{Sensitivity of the test, so the values on the \eqn{y}-axis of
#' a ROC plot.}
#' \item{ppvs}{Positive predictive values for each cutoff.}
#' \item{npvs}{Negative predictive values for each cutoff.}
#' \item{cutoffs}{Cutoffs used (basically the pooled marker values of
#' cases and controls).}
#' \item{res.mat}{Collects the above quantities in a matrix, including
#' Wilson confidence intervals, computed at at confidence level \code{conf.level}.}
#' \item{auc}{Area under the ROC curve. This is equal to the value of the
#' Mann-Whitney test statistic.}
#' \item{auc.var}{Variance of AUC.}
#' \item{auc.var.norm}{Variance of AUC if data is assumed to come from a
#' bivariate normal distribution.}
#' \item{lowCI}{Lower limit of Wald confidence interval.}
#' \item{upCI}{Upper limit of Wald confidence interval.}
#' \item{logitLowCI}{Lower limit of a Wald confidence interval received on
#' logit scale.}
#' \item{logitUpCI}{Upper limit of a Wald confidence interval received on
#' logit scale.}
#' @note The confidence intervals are only valid if observations are
#' \emph{independent}.
#' @author Kaspar Rufibach \email{kaspar.rufibach@@gmail.com} and Andrea Riebler.
#' @seealso \code{\link{confIntAUCbinorm}}. Similar functionality is provided in the package \pkg{ROCR}.
#' @references The original reference for the computation of the confidence interval is:
#' 
#' Hanley, J.A. and McNeil, B.J. (1982). The meaning and use of the area under
#' the curve. \emph{Radiology}, \bold{143}, 29--36.
#' 
#' See also:
#' 
#' Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for
#' classification and prediction}. Oxford University Press.
#' @keywords htest
#' @examples
#'
#' ## simulate data
#' ## --------------
#' set.seed(1977)
#' controls <- rnorm(n = 50)
#' cases <-  rnorm(n = 40, mean = 0.5, sd = 1.5)
#' 
#' ## summary of ROC curve
#' ## --------------
#' res <- summaryROC(cases, controls, conf.level = 0.95)
#' 
#' 
#' ## alternative bootstrap CI for AUC
#' ## --------------
#' resBinormBoot <- confIntAUCbinorm(cases = cases, controls = controls, 
#'                                   conf.level = 0.95, replicates = 1000,
#'                                   grid = 100)
#' 
#' ## alternative bootstrap CI for AUC 
#' ## --------------
#' resBinormWald <- confIntAUCbinorm(cases, controls, conf.level = 0.95,
#'                                   method = "wald")
#' 
#' ## display results
#' ## --------------
#' str(res)
#' str(resBinormBoot)
#' str(resBinormWald)
#' 
#' ## plot ROC curve
#' ## --------------
#' plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
#'      xlab = "1 - specificity", ylab = "sensitivity", pty = 's')
#' segments(0, 0, 1, 1, lty = 2)
#' lines(res$x.val, res$y.val, type = 'l', col = 2, lwd = 2, lty = 2)
#' lines(resBinormBoot$x.val, resBinormBoot$y.val, type = 'l',
#'       col = 4, lwd = 2, lty = 2)
#' 
#' @export
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

    stopifnot(is.numeric(cases), length(cases) >= 1, is.finite(cases),
              is.numeric(controls), length(controls) >= 1, is.finite(controls),
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level), 0 < conf.level, conf.level < 1)
    
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
    npvs.ci <- ppvs.ci <- sens.ci <- spec.ci <- matrix(NA, nrow = n, ncol = 2)
    
    for (i in 1:n){
        spec.ci[i, ] <- wilson(x = tns[i], n = (tns + fps)[i],
                               conf.level = conf.level)[c(1, 3)]
        sens.ci[i, ] <- wilson(x = tps[i], n = (tps + fns)[i],
                               conf.level = conf.level)[c(1, 3)]
        ppvs.ci[i, ] <- wilson(x = tps[i], n = (tps + fps)[i],
                               conf.level = conf.level)[c(1, 3)]
        npvs.ci[i, ] <- wilson(x = tns[i], n = (fns + tns)[i],
                               conf.level = conf.level)[c(1, 3)]
    }
    
    # compute q2 and q1
    q2 <- sum(n3tmp * (tns ^ 2 + tns * n1tmp + 1/3 * n1tmp ^ 2)) /
        (n2 * n1 ^ 2)
    q1 <- sum(n1tmp *(tps ^ 2 + tps * n3tmp + 1/3 * n3tmp ^ 2)) /
        (n1 * n2 ^ 2)

    # estimate auc according to Hanley and McNeil (1982)
    auc <- sum(n1tmp * tps + 0.5 * n1tmp * n3tmp) / (n1 * n2)

    # estimate AUC as test statistic of Wilcoxon test
    auc <- as.numeric(wilcox.test(cases, controls, exact = FALSE)$statistic /
                               (n1 * n2))
    auc.var <- (n1 * n2) ^ (-1) * (auc * (1 - auc) + (n2 - 1) * (q1 - auc ^ 2) +
                                   (n1 - 1) * (q2 - auc ^ 2))

    # compute variance of AUC when assuming a normal model for measurements
    q1norm <- auc / (2 - auc)
    q2norm <- 2 * auc ^ 2 / (1 + auc)
    auc.var.norm <- (n1 * n2) ^ (-1) * (auc * (1 - auc) + (n2 - 1) * (q1norm - auc ^ 2) +
                                        (n1 - 1) * (q2norm - auc ^ 2))

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
    res.mat <- cbind(cutoffs, 1 - x.val, spec.ci, y.val, sens.ci,
                     ppvs, ppvs.ci, npvs, npvs.ci)
    colnames(res.mat) <- c("cutoff", "specificity", "CIspeclow", "CIspecup",
                           "sensitivity", "CIsenslow", "CIsensup", "PPV", "CIPPVlow",
                           "CIPPVup", "NPV", "CINPVlow", "CINPVup")

    list("x.val" = x.val, "y.val" = y.val, "ppvs" = ppvs, "npvs" = npvs,
         "cutoffs" = cutoffs, "res.mat" = res.mat, "auc" = auc, 
         "auc.var" = auc.var, "auc.var.norm" = auc.var.norm, "lowCI" = lowCI,
         "upCI" = upCI, "logitLowCI" = logLowCI, "logitUpCI" = logUpCI)
}
