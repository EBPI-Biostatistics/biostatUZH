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
#' @param status Numeric vector of censoring indicators, 1 for event, 0 for censored.
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
#' time <- leukemia[, 1]; status <- leukemia[, 2]; group <- leukemia[, 3]
#' survReport(time = time, status = status)
#' survReport(time = time, status = status, group = group)
#' survReport(time = time, status = status, group = group, stmt.placement = "bottomright",
#'            output = "plain", main = "Acute Myeloid Leukemia")
#' survReport(time = time, status = status, group = group, stmt.placement = "bottomright",
#'            output = "text", main = "AML")
#' survReport(time = time, status = status, group = group, stmt.placement = "subtitle",
#'            output = "text", main = "Acute Myeloid Leukemia")
#' 
#' ## use Larynx data from the 'survival' package
#' data(larynx)
#' larynx$stage <- factor(x = larynx$stage, 1:4,  c("I", "II", "III", "IV"))
#' with(larynx, survReport(time = time, status = death, group = stage, output = "text",
#'                         stmt.placement = "bottomright", main = "Larynx"))
#' with(larynx, survReport(time = time, status = death, group = stage, output = "text",
#'                         stmt.placement = "subtitle", main = ""))
#' with(larynx, survReport(time = time, status = death, group = stage, output = "plain",
#'                         stmt.placement = "topright", legend.placement = "bottomright"))
#'
#' @import survival
#' @export
survReport <- function(time, status, group = NULL, 
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
              is.numeric(status),
              length(status) == length(time))
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
        fm <- survival::Surv(time, status) ~ 1
        np <- 1
        medians <- quantileKM(time = time, event = status, group = NA, conf.level = conf.level)$quantities[3:5]
        med.stmt <- paste(round(medians[1], digits), " (", round(medians[2], digits), " - ",
                          round(medians[3], digits), ")", sep = "")
    } else {
        fm <- survival::Surv(time, status) ~ group
        np <- length(table(group))
        medians <- quantileKM(time = time, event = status, group = group, conf.level = conf.level)$quantities[,3:5]
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
