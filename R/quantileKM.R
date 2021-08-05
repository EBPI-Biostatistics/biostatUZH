    # Compute basic quantities for Kaplan-Meier estimate, using
    # functions provided by library(survival)
    #
    # Input: 
    #   - time:    event times
    #   - event:   censoring indicator (0 = censored, 1 = event)
    #   - group:   grouping factor, if NA --> only one group
    #   - quant:   quantile to be computed, defaults to median
    #   - alpha:   confidence level for confidence intervals
    #
    # Kaspar Rufibach, October 2008
    # Update Leonhard Held, November 2016 




#' Compute arbitrary quantile of a Kaplan-Meier estimate of a survival function
#' 
#' In survival analysis, often some quantiles (mainly the median) of an
#' estimated survival function are of interest. The survfit function in the
#' 'survival' packages of version \eqn{< 2.35} computed median time to event.
#' However, the output of that function was such that the median is not
#' accessible by the user (it can only be read off the output). This function
#' makes the median and any other quantile accessible by the user.
#' 
#' 
#' @param time Event times, censored or observed.
#' @param event Censoring indicator, 1 for event, 0 for censored.
#' @param group Quantiles can be computed for several survival curves, defined
#' by \code{group}.
#' @param quant Quantile to be computed. Real number in \eqn{[0, 1]}.
#' @param conf.level Significance level for confidence interval for the time to
#' event quantile.
#' @param conftype Type of confidence interval to be computed. For possible
#' choices see above, and for specifications regarding the different options
#' confer the help file of the function \code{survfit}.
#' @param conflower Controls modified lower limits to the curve, the upper
#' limit remains unchanged. See \code{survfit} for details.
#' @return \item{n}{Number of observations used.} \item{events}{Number of
#' events.} \item{quantile}{Quantile estimate.} \item{lower.ci}{Lower limit of
#' confidence interval.} \item{upper.ci}{Upper limit of confidence interval.}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @seealso Partly based on the function \code{survfit}.
#' @references Computation of confidence intervals is done according to
#' 
#' Brookmeyer, R. and Crowley, J. (1982).  A Confidence Interval for the Median
#' Survival Time.  \emph{Biometrics}, \bold{38}, 29--41.
#' @keywords htest survival
#' @examples
#' 
#' ## use Acute Myelogenous Leukemia survival data contained in package 'survival'
#' time <- leukemia[, 1]; status <- leukemia[, 2]; x <- as.factor(leukemia[, 3])
#' plot(survfit(Surv(time, status) ~ x, conf.type = "none"), mark = "/", col = 1:2)
#' 
#' ## median time to event
#' qKM <- quantileKM(time, status, group = x, quant = 0.5, conf.level = 0.95, 
#'     conftype = "log")
#' qKM
#' 
#' ## extract results
#' qKM$quantities
#' 
#' ## comparison to standard function (median time to event not accessible by user)
#' quant.surv <- survfit(Surv(time, status) ~ x, conf.int = 0.95)
#' quant.surv
#' 
#' ## compute 0.25 quantile
#' qKM2 <- quantileKM(time, status, group = x, quant = 0.25, conf.level = 0.95, 
#'     conftype = "log")
#' qKM2
#'
#' @import survival
#' @export
quantileKM <- function(time, event, group = NA, quant = 0.5, conf.level = 0.95, conftype = c("log","log-log","plain","none")[2],
                       conflower = c("usual", "peto", "modified")[1]){

    
    ## initialization
    alpha <- 1 - conf.level
    s.obj <- Surv(time, event)
    
    if (identical(group, NA) == FALSE){
        group.f <- as.factor(group)
        n.level <- length(levels(group.f))
        quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
        sdiff <- survdiff(s.obj ~ group.f)
        
        # p-value log-rank test
        p.val <- 1 - pchisq(sdiff$chisq, df = n.level - 1)
        quant.surv <- survfit(s.obj ~ group.f, conf.int = 1 - alpha, conf.type = conftype, conf.lower = conflower)

        ## quantile of event time, incl. confidence intervals
        for (j in 1:n.level){
            tmp <- summary(quant.surv[j])
            quant.mat[j, 1] <- quant.surv[j]$n
            quant.mat[j, 2] <- sum(quant.surv[j]$n.event)
        
            ## add Inf to omit warnings in case quantile is not reached by 
            ## survival curve (or pointwise ci curve)
            quant.mat[j, 3] <- min(Inf, tmp$time[tmp$surv <= quant], na.rm = TRUE)
            quant.mat[j, 4] <- min(Inf, tmp$time[tmp$lower <= quant], na.rm = TRUE)
            quant.mat[j, 5] <- min(Inf, tmp$time[tmp$upper <= quant], na.rm = TRUE)
        }        
        
    dimnames(quant.mat)[[1]] <- paste("group = ", levels(group.f), sep = "")
    }
        
    if (identical(group, NA) == TRUE){
        p.val <- NA
        n.level <- 1
        quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
        quant.surv <- survfit(s.obj ~ 1, conf.int = 1 - alpha, conf.type = conftype, conf.lower = conflower)
        tmp <- quant.surv
        quant.mat[1, 1] <- tmp$n
        quant.mat[1, 2] <- sum(tmp$n.event)
        
        ## add Inf to omit warnings in case quantile is not reached by 
        ## survival curve
        quant.mat[1, 3] <- min(Inf, tmp$time[tmp$surv <= quant], na.rm = TRUE)
        quant.mat[1, 4] <- min(Inf, tmp$time[tmp$lower <= quant], na.rm = TRUE)
        quant.mat[1, 5] <- min(Inf, tmp$time[tmp$upper <= quant], na.rm = TRUE)        
        }    
        
dimnames(quant.mat)[[2]] <- c("n", "events", "quantile", "lower.ci", "upper.ci")      
return(list("quantities" = quant.mat, "p.val" = p.val))
}
