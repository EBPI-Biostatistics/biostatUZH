#' Survival Analysis Results for Reports
#' 
#' The analysis of survival data most often requires the production of
#' Kaplan-Meier plots, estimated median survival time, and univariate estimates
#' of hazard ratios (with associated confidence intervals and p-values). Often,
#' this is performed separately for several possible risk factors. This
#' function simplifies these tasks by providing a unified interface to
#' simultaneously produce graphics, estimate median survival time, and compute
#' hazard ratios with associated statistics from the Cox proportional hazards
#' model.
#' 
#' 
#' @param time Numeric vector of event times, censored or observed.
#' @param event Numeric vector of censoring indicators, 1 for event, 0 for censored.
#' @param group Kaplan-Meier plots, median survival and hazard rates will be
#' computed for each \code{group}. If \code{group} is missing, the Kaplan-Meier plots
#' and median survival will be produced for the entire dataset.
#' @param stmt.placement Placement of the statement giving hazard ratios, confidence
#' intervals and p-values in the plot. Ignored if \code{group = NA}.
#' @param legend.placement Placement of the legend in the plot. Ignored if
#' \code{group = NA}. A warning will be produced if \code{stmt.placement} equals \code{legend.placemet}.
#' @param output Type of output, "plain" gives plaintext output for median
#' survival and hazard ratios suitable for the creation of tables etc., while
#' "text" gives output suitable for the text of reports. "text"-type output
#' is always printed in the plots unless \code{stmt.placement} is "none".
#' @param main Title of the plot. 
#' @param xlab Label for the x-axis of the plot. 
#' @param ylab Label for the y-axis of the plot. 
#' @param labels Labels used to produce the legend and the contrasts
#' printed in the hazard ratio statements.
#' @param digits Number of significant digits used in rounding for text-type
#' output.
#' @param conf.level Significance level for confidence interval for the median
#' survival and hazard ratio estimates.
#' @return The Kaplan-Meier curve is plotted and a list with the following elements is returned:
#' \item{med}{Matrix of median survival, and corresponding confidence intervals.
#' If \code{output = "text"}, this is a vector of text giving
#' median survival estimates and confidence intervals.}
#' \item{events}{Matrix of hazard ratios, corresponding confidence intervals, and p-values from the
#' Wald test. If \code{output = "text"}, this is a character vector.}
#' @author Sarah R. Haile
#' @seealso \code{survfit}, \code{coxph}, \code{quantileKM}
#' @keywords survival
#' @examples
#' 
#' ## use Acute Myelogenous Leukemia survival data from the 'survival' package
#' time <- leukemia[, 1]; event <- leukemia[, 2]; group <- leukemia[, 3]
#' survReport(time = time, event = event)
#' survReport(time = time, event = event, group = group)
#' survReport(time = time, event = event, group = group, stmt.placement = "bottomright",
#'            output = "plain", main = "Acute Myeloid Leukemia")
#' survReport(time = time, event = event, group = group, stmt.placement = "bottomright",
#'            output = "text", main = "AML")
#' survReport(time = time, event = event, group = group, stmt.placement = "subtitle",
#'            output = "text", main = "Acute Myeloid Leukemia")
#' 
#' ## use Larynx data from the 'survival' package
#' data(larynx)
#' larynx$stage <- factor(x = larynx$stage, 1:4,  c("I", "II", "III", "IV"))
#' with(larynx, survReport(time = time, event = death, group = stage, output = "text",
#'                         stmt.placement = "bottomright", main = "Larynx"))
#' with(larynx, survReport(time = time, event = death, group = stage, output = "text",
#'                         stmt.placement = "subtitle", main = ""))
#' with(larynx, survReport(time = time, event = death, group = stage, output = "plain",
#'                         stmt.placement = "topright", legend.placement = "bottomright"))
#'
#' @import survival
#' @export
survReport <- function(time, event, group = NULL, 
                       stmt.placement = c("bottomleft", "bottomright", "topright", "subtitle", "none"), 
                       legend.placement = c("topright", "bottomleft", "bottomright", "none"),
                       output = c("text", "plain"),
                       main = "",
                       xlab = "Time", 
                       ylab = "Survival",
                       labels = levels(group),
                       digits = 2,
                       conf.level = 0.95){

    stopifnot(is.numeric(time),
              length(time) > 0,
              is.numeric(event),
              length(event) == length(time))
    if(!is.null(group)){
        group <- as.factor(group)
        stopifnot(is.finite(group),
                  length(group) == length(time))
    }
    stopifnot(!is.null(stmt.placement))
    stmt.placement <- match.arg(stmt.placement)
    stopifnot(!is.null(legend.placement))
    legend.placement <- match.arg(legend.placement)
    stopifnot(!is.null(output))
    output <- match.arg(output)
    stopifnot(is.character(main),
              length(main) == 1,
              is.character(xlab),
              length(xlab) == 1,
              is.character(ylab),
              length(ylab) == 1)
    if(!is.null(group)){
        stopifnot(is.character(labels),
                  length(labels) == length(unique(group)))
    }
    stopifnot(is.numeric(digits),
              length(digits) == 1,
              is.finite(digits),
              digits >= 0,
              is.numeric(conf.level),
              length(conf.level) == 1,
              is.finite(conf.level),
              0 <= conf.level, conf.level <= 1)
    
    if(stmt.placement == legend.placement && stmt.placement!="none")
        warning("Legend and Hazard Ratio Statement will be overplotted! Change legend or statment placement option.")

    if (is.null(group)) {
        fm <- survival::Surv(time, event) ~ 1
        np <- 1
        medians <- quantileKM(time = time, event = event, group = NA, conf.level = conf.level)$quantities[3:5]
        med.stmt <- paste(round(medians[1], digits), " (", round(medians[2], digits), " - ",
                          round(medians[3], digits), ")", sep = "")
    } else {
        fm <- survival::Surv(time, event) ~ group
        np <- length(table(group))
        medians <- quantileKM(time = time, event = event, group = group, conf.level = conf.level)$quantities[,3:5]
        med.stmt <- paste(round(medians[,1], digits), " (", round(medians[,2], digits), " - ",
                          round(medians[,3], digits), ")", sep = "")
    }
    
    sv <- survival::survfit(formula = fm)
    tab <- summary(object = sv)$table
    
    plot(sv, col = 1:np, lty = 1:np, conf.int = FALSE, main = main, xlab = xlab, ylab = ylab)
    if(np > 1){
        cph <- coxph(fm)
        cph.var <- cph$var
        if(np > 2)
            cph.var <- diag(cph.var)
        z <- coef(cph) / sqrt(cph.var)
        if(np == 2){
            pval <- 1 - pchisq(survival::survdiff(fm)$chisq,
                               df = length(survival::survdiff(fm)$n) - 1)
        } else {
            pval <- pnorm(-abs(z)) * 2
        }
        contra <- paste(labels[1], " vs. ", labels[2:np], sep = "")
        hr.txt <- cbind(exp(coef(cph)), exp(confint(cph, level = conf.level)), pval)
        rownames(hr.txt) <- contra
        
        hr.stmt <- paste(contra, ": HR ", round(hr.txt[,1], digits), " (", round(hr.txt[,2], digits),
                         " - ", round(hr.txt[,3], digits), "), p ", ifelse(hr.txt[,4] < 0.0001, "", "= "),
                         format.pval(hr.txt[,4], digits = digits, eps = 0.0001), sep = "") 
        hr.stmt2 <- hr.stmt
        if(!stmt.placement %in% c("topright"))
            hr.stmt2 <- rev(hr.stmt)
        if(!stmt.placement %in% c("subtitle", "none")){
            mult <- switch(stmt.placement, bottomleft = 1 / 50, bottomright = 1 - 1 / 50, topright = 1 - 1/50)
            algn <- switch(stmt.placement, bottomleft = 4, bottomright = 2, topright = 2)
            incr <- (1:(np - 1)) / 20
            vert <- switch(stmt.placement, bottomleft = 0 + incr, bottomright = 0 + incr, topright = 1 - incr, none = 0)
            text(x = max(sv$time) * mult, y = vert, labels = hr.stmt2, pos = algn)
        }
        if(stmt.placement == "subtitle"){    
            mtext(hr.stmt2, side = 3, line = 0:(np - 2))
        }
        if(legend.placement != "none")
            legend(legend.placement, labels, col = 1:np, lty = 1:np, bty = "n")
    }
    
    if(output == "plain"){
        tmp.out <- list(med = medians)
        if(np > 1)
            tmp.out$hr <- hr.txt
    } else {
        tmp.out <- list(med = med.stmt)
        if(np > 1)
            tmp.out$hr <- hr.stmt
    }
    tmp.out
}



#' Quantile of a Kaplan-Meier estimate of a survival function
#' 
#' In survival analysis, often quantiles (e.g., the median) of an
#' estimated survival function are of interest. The \code{\link[survival]{survfit}}
#' function in the  'survival' packages of version \eqn{< 2.35} computes median time to event.
#' However, the output of that function is such that the median is not
#' conveniently accessible by the user. This function
#' makes the median and any other quantile accessible by the user.
#' 
#' @param time Event times, censored or observed.
#' @param event Censoring indicator, 1 for event, 0 for censored.
#' @param group Indicates groups to compute individual quantile for each group.
#' Default is \code{NULL}.
#' @param quantile Quantile to be computed. Real number in \eqn{[0, 1]}.
#' @param conf.level Significance level for confidence interval for the time to
#' event quantile.
#' @param conf.type Type of confidence interval to be computed. For possible
#' choices see above, and for specifications regarding the different options
#' see \code{\link[survival]{survfit}}.
#' @return
#' A list containing a matrix with columns:
#' \item{n}{Number of observations used.}
#' \item{events}{Number of events.}
#' \item{quantile}{Quantile estimate.}
#' \item{lower.ci}{Lower limit of confidence interval.}
#' \item{upper.ci}{Upper limit of confidence interval.}\cr
#' and the p-value of the test for the difference between the survival curves;
#' see \code{\link[survival]{survdiff}} for more information.
#' @author Kaspar Rufibach
#' @note The function is based on the function \code{\link[survival]{survfit}}.
#' @seealso \code{\link[survival]{survfit}}, \code{\link[survival]{survdiff}}
#' @references Confidence intervals are computed according to
#' 
#' Brookmeyer, R. and Crowley, J. (1982).  A Confidence Interval for the Median
#' Survival Time.  \emph{Biometrics}, \bold{38}, 29--41.
#' @keywords htest survival
#' @examples
#' 
#' ## use Acute Myelogenous Leukemia survival data contained in package 'survival'
#' time <- leukemia[, 1]
#' event <- leukemia[, 2]
#' group <- as.factor(leukemia[, 3])
#' plot(survfit(Surv(time, event) ~ group, conf.type = "none"), mark = "/", col = 1:2)
#' 
#' ## median time to event
#' quantileKM(time, event, group = group, quantile = 0.5, conf.level = 0.95, 
#'            conf.type = "log")
#'  
#' ## comparison to standard function (median time to event not accessible by user)
#' survfit(Surv(time, event) ~ group, conf.int = 0.95)
#' 
#' ## compute 0.25 quantile
#' quantileKM(time = time, event = event, group = group, quantile = 0.25,
#'            conf.type = "log")
#'
#' @import survival
#' @export
quantileKM <- function(time, event, group = NULL, quantile = 0.5, conf.level = 0.95,
                       conf.type = c("log-log", "log", "plain","none")){
    
    stopifnot(is.numeric(time),
              length(time) > 0,
              is.numeric(event),
              length(event) == length(time))
    if(!is.null(group)){
        group <- as.factor(group)
        stopifnot(is.finite(group),
                  length(group) == length(time))
    }
    stopifnot(is.numeric(quantile), length(quantile) == 1,
              is.finite(quantile),
              0 <= quantile, quantile <= 1,
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    stopifnot(!is.null(conf.type))
    conf.type <- match.arg(conf.type)
    
    s.obj <- Surv(time, event)
    
    if (!is.null(group)) {
        group.f <- as.factor(group)
        n.level <- length(levels(group.f))
        quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
        sdiff <- survdiff(s.obj ~ group.f)
        
        ## p-value log-rank test
        p.val <- 1 - pchisq(sdiff$chisq, df = n.level - 1)
        quant.surv <- survfit(s.obj ~ group.f, conf.int = conf.level,
                              conf.type = conf.type)
        
        ## quantile of event time, incl. confidence intervals
        for (j in 1:n.level){
            tmp <- summary(quant.surv[j])
            quant.mat[j, 1] <- quant.surv[j]$n
            quant.mat[j, 2] <- sum(quant.surv[j]$n.event)
            
            ## add Inf to omit warnings in case quantile is not reached by 
            ## survival curve (or pointwise ci curve)
            quant.mat[j, 3] <- min(Inf, tmp$time[tmp$surv <= quantile], na.rm = TRUE)
            quant.mat[j, 4] <- min(Inf, tmp$time[tmp$lower <= quantile], na.rm = TRUE)
            quant.mat[j, 5] <- min(Inf, tmp$time[tmp$upper <= quantile], na.rm = TRUE)
        }        
        
        dimnames(quant.mat)[[1]] <- paste("group = ", levels(group.f), sep = "")
    }
    
    if (is.null(group)) {
        p.val <- NA
        n.level <- 1
        quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
        quant.surv <- survfit(s.obj ~ 1, conf.int = conf.level, conf.type = conf.type)
        tmp <- quant.surv
        quant.mat[1, 1] <- tmp$n
        quant.mat[1, 2] <- sum(tmp$n.event)
        
        ## add Inf to omit warnings in case quantile is not reached by 
        ## survival curve
        quant.mat[1, 3] <- min(Inf, tmp$time[tmp$surv <= quantile], na.rm = TRUE)
        quant.mat[1, 4] <- min(Inf, tmp$time[tmp$lower <= quantile], na.rm = TRUE)
        quant.mat[1, 5] <- min(Inf, tmp$time[tmp$upper <= quantile], na.rm = TRUE)        
    }    
    
    dimnames(quant.mat)[[2]] <- c("n", "events", "quantile", "lower.ci", "upper.ci")      
    list("quantities" = quant.mat, "p.val" = p.val)
}
