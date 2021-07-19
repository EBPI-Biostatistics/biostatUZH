#' Function to compute ROC curve and (bootstrap) confidence interval for AUC in
#' the binormal model
#' 
#' This function computes for values of a continuous variable of a group of
#' cases and a group of controls the ROC curve in the binormal model.
#' Additionally, the AUC including a confidence interval is provided, with
#' bootstrap or with Wald, assuming the variances are equal (the assumption of
#' unequal variances has yet to be implemented).
#' 
#' 
#' @aliases confIntAUCbinorm AUCbinorm
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @param ci.method Method of calculating the confidence interval. Implemented
#' are \code{boot} for bootstrap and \code{wald} for the Wald based confidence
#' interval.
#' @param conf.level Confidence level for confidence interval.
#' @param replicates Number of boostrap replicates. Only needed if
#' \code{ci.method = 'boot'}.
#' @param grid The interval \eqn{[0, 1]} is split in \code{grid} intervals to
#' compute the ROC curve on. Only needed if \code{ci.method = 'boot'}.
#' @param var.equal Logical, are the variances assumed to be equal or not (only
#' \code{TRUE} implemented). Only needed if \code{ci.method = 'wald'}.
#' @return The results for \code{ci.method = 'boot'} are \item{a}{The values of
#' a.} \item{b}{The values of b.} \item{x.val}{The values on the \eqn{x}-axis
#' of a ROC plot.} \item{y.val}{The values on the \eqn{y}-axis of a ROC plot,
#' i.e. the binormal ROC curve.} \item{auc}{Area under the ROC curve. }
#' \item{lowBootCI}{Lower limit of bootstrap confidence interval.}
#' \item{upBootCI}{Upper limit of bootstrap confidence interval.}
#' 
#' The results for \code{ci.method = 'wald'} are \item{a}{The values of a.}
#' \item{b}{The values of b.} \item{auc}{Area under the ROC curve.}
#' \item{lowCI}{Lower limit of confidence interval.} \item{upCI}{Upper limit of
#' confidence interval.}
#' @note The Wald approach is documented in a vignette: see
#' \code{vignette("aucbinormal")}
#' @author Kaspar Rufibach (\code{ci.method = 'boot'}) and Leonhard Held
#' (\code{ci.method = 'wald'})
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford: Oxford University Press.
#' @keywords htest
#' @examples
#' 
#' set.seed(1977)
#' ns <- c(50, 40)
#' truth <- c(rep(0, ns[1]), rep(1, ns[2]))
#' estimates <- c(rnorm(ns[1]), rnorm(ns[2], mean = 0.5, sd = 1.5))
#' cases <- estimates[truth == 1]
#' controls <- estimates[truth == 0]
#' res <- summaryROC(cases, controls, conf.level = 0.95)
#' 
#' ## ci.method = 'boot'
#' ## --------------
#' resBinorm <- confIntAUCbinorm(cases, controls, 
#'     conf.level = 0.95, replicates = 1000, grid = 100)
#' 
#' # display results
#' res
#' resBinorm
#' 
#' ## plot ROC curve
#' ## --------------
#' plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
#'     xlab = "1 - specificity", ylab = "sensitivity", pty = 's')
#' segments(0, 0, 1, 1, lty = 2)
#' lines(res$x.val, res$y.val, type = 'l', col = 2, lwd = 2, lty = 2)
#' lines(resBinorm$x.val, resBinorm$y.val, type = 'l', col = 4, lwd = 2, lty = 2)
#' 
#' ## ci.method = 'wald'
#' ## --------------
#' resBinorm <- confIntAUCbinorm(cases, controls, conf.level = 0.95, ci.method = "wald")
#' 
confIntAUCbinorm <- function(cases, controls, conf.level = 0.95, replicates =
                             1000, grid = 100, ci.method = c("boot", "wald"),
                             var.equal = TRUE){ 

    ## decide CI method
    ci.method <- match.arg(ci.method)
    
    alpha <- 1 - conf.level

    ## CI with bootstrap
    ## -------------------
    if (ci.method == "boot")
    {
        ## compute bootstrap CI for binormal AUC
        ## data in two columns: y (measurements), status (0 = control, 1 = case)
        
        ts <- seq(0, 1, by = 1 / grid)		
        a <- (mean(cases) - mean(controls)) / sd(cases)
        b <- sd(controls) / sd(cases)	
        biroc <- pnorm(a + b * qnorm(ts))

        dat <- cbind(c(controls, cases), c(rep(0, length(controls)), rep(1, length(cases))))

        ## compute bootstrap ci
        res <- boot(dat, AUCbinorm, R = replicates)
        aucbin <- boot.ci(res,type = "bca")

        return(list("a" = a, "b" = b, "x.val" = ts, "y.val" = biroc,
                    "auc" = res$t0, "lowBootCI" = aucbin$bca[4],
                    "upBootCI" = aucbin$bca[5]))  
    }

    ## CI 
    ## -------------------
    if(ci.method == "wald")
    {
        n0 <- length(controls)
        n1 <- length(cases)
        
        mu0 <- mean(controls)
        mu1 <- mean(cases)
        
        if(var.equal)
        {
            s0 <- s1 <- sd(c(cases, controls))
            a.se <- sqrt(1/n0 + 1/n1)
            auc.se <- a.se
        } else {
            stop("var.equal = FALSE not implemented for ci.method = 'wald'")
            
            ## s0 <- sd(controls)
            ## s1 <- sd(cases)
        }

        a <- (mu1 - mu0)/s1
        b <- s0/s1
        auc <- pnorm(a/(sqrt(1 + b^2)))
        lowCI <- pnorm((a - qnorm(1-alpha/2) * auc.se) / (sqrt(1 + b^2)))
        upCI <- pnorm((a + qnorm(1-alpha/2) * auc.se) / (sqrt(1 + b^2)))
        
        return(list("a" = a, "b" = b, "auc" = auc, "lowCI" = lowCI, "upCI" = upCI))
    }


}
