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
#' @param ftime Event times, censored or observed.
#' @param fstatus Censoring indicator, 1 for event, 0 for censored.
#' @param fgroup Kaplan-Meier plots, median survival and hazard rates will be
#' computed by \code{fgroup}. If \code{fgroup} is missing, Kaplan-Meier plots
#' and median survival will be produced for the entire dataset.
#' @param stmt.pl Placement of the statement giving hazard ratios, confidence
#' intervals and p-values in the plot. Ignored if \code{fgroup = NA}.
#' @param legend.pl Placement of the legend in the plot. Ignored if
#' \code{fgroup = NA}. A warning will be produced if \code{stmt.pl =
#' legend.pl}.
#' @param output Type of output, ``plain'' gives plaintext output for median
#' survival and hazard ratios suitable for the creation of tables etc., while
#' ``text'' gives output suitable for the text of reports. ``text''-type output
#' is always printed in the plots unless \code{stmt.pl="none"}.
#' @param maintitle Main title of the plot. Corresponds to the \code{main}
#' option of \code{plot.survfit}.
#' @param ylbl Label for the y-axis of the plot. Corresponds to the \code{ylab}
#' option of \code{plot.survfit}.
#' @param xlbl Label for the x-axis of the plot. Corresponds to the \code{xlab}
#' option of \code{plot.survfit}.
#' @param lbls Labels used to produce the legend, as well as the contrasts
#' printed in the hazard ratio statements.
#' @param dig Number of significant digits used in rounding for ``text''-type
#' output.
#' @param conf.level Significance level for confidence interval for the median
#' survival and hazard ratio estimates.
#' @return The Kaplan-Meier plot is printed, as well as a list containing:
#' \item{med}{Matrix of median survival, and corresponding confidence
#' intervals. If \code{output = "text"}, this will be a vector of text giving
#' median survival estimates and confidence intervals.} \item{events}{Matrix of
#' hazard ratios, corresponding confidence intervals, and p-values from the
#' Wald test. If \code{output = "text"}, this will be a vector of text.}
#' @author Sarah R. Haile
#' @seealso \code{survfit}, \code{coxph}, \code{quantileKM}
#' @keywords survival
#' @examples
#' 
#' ## use Acute Myelogenous Leukemia survival data contained in package 'survival'
#' time <- leukemia[, 1]; status <- leukemia[, 2]; x <- leukemia[, 3]
#' survrep(time, status, NA)
#' survrep(time, status, x)
#' survrep(time, status, x, stmt.pl = "bottomright", out = "plain", main = "Acute Myeloid Leukemia")
#' survrep(time, status, x, stmt.pl = "bottomright", out = "text", main = "AML")
#' survrep(time, status, x, stmt.pl = "subtitle", out = "text", main = "Acute Myeloid Leukemia")
#' 
#' ## using Larynx data contained in package 'survival'
#' ## to show behavior with risk factors with >2 levels
#' data(larynx)
#' larynx$stage <- factor(larynx$stage, 1:4,  c("I", "II", "III", "IV"))
#' with(larynx, survrep(time, death, stage, out = "t",
#'                      stmt.pl = "bottomright", main = "Larynx"))
#' with(larynx, survrep(time, death, stage, out = "t",
#'                      stmt.pl = "subtitle", main = ""))
#' with(larynx, survrep(time, death, stage, out = "p",
#'                      stmt.pl = "topright", legend.pl = "bottomright"))
#' 
#' @export
survrep <- function(ftime, fstatus, fgroup = NA, 
        stmt.pl = c("bottomleft", "bottomright", "topright", "subtitle", "none")[1], 
        legend.pl = c("bottomleft", "bottomright", "topright", "none")[3],
        output = c("plain", "text")[2],
        maintitle = "",
        ylbl = "Survival",
        xlbl = "Time", 
        lbls = levels(fgroup),
        dig = 2, conf.level = 0.95){

if (identical(fgroup, NA) == TRUE) {
    fm <- Surv(ftime, fstatus) ~ 1
    np <- 1
    
} else {
    fm <- Surv(ftime, fstatus) ~ fgroup
    np <- length(table(fgroup))
}

out <- match.arg(output, c("plain", "text"))
stmt.placement <- match.arg(stmt.pl, c("bottomleft", "bottomright", "topright", "subtitle",  "none"))
legend.placement <- match.arg(legend.pl, c("bottomleft", "bottomright", "topright", "none"))
if((stmt.placement == legend.placement) & stmt.placement!="none") warning("Legend and Hazard Ratio Statement will be overplotted! Change legend or statment placement option.")


sv <- survfit(fm)
tab <- summary(sv)$table

if(np==1){
medians <- quantileKM(ftime, fstatus, fgroup, conf.level = conf.level)$quantities[3:5]
med.stmt <- paste(round(medians[1], dig), " (", round(medians[2], dig), " - ", round(medians[3], dig), ")", sep = "")
} else {
medians <- quantileKM(ftime, fstatus, fgroup, conf.level = conf.level)$quantities[,3:5]
med.stmt <- paste(round(medians[,1], dig), " (", round(medians[,2], dig), " - ", round(medians[,3], dig), ")", sep = "")
}

plot(sv, col = 1:np, lty = 1:np, conf.int = FALSE, 
        main = maintitle, xlab = xlbl, ylab = ylbl)
summary(sv)
if(np>1){
    cph <- coxph(fm)
    cph.var <- cph$var
    if(np>2) cph.var <- diag(cph.var)
    z <- coef(cph)/sqrt(cph.var)
    if(np==2){
            pval <- 1 - pchisq(survdiff(fm)$chisq, df = length(survdiff(fm)$n) - 1)
    } else  pval <- pnorm(-abs(z))*2
    contra <- paste(lbls[1], " vs. ", lbls[2:np], sep = "")
    hr.txt <- cbind(exp(coef(cph)), exp(confint(cph, level = conf.level)), pval)
    rownames(hr.txt) <- contra

    hr.stmt <- paste(contra, ": HR ", round(hr.txt[,1], dig), " (", round(hr.txt[,2], dig), " - ", round(hr.txt[,3], dig), "), p ", ifelse(hr.txt[,4]<0.0001, "", "= "), format.pval(hr.txt[,4], digits = dig, eps = 0.0001), sep = "") 
    hr.stmt2 <- hr.stmt
    if(!stmt.placement %in% c("topright")) hr.stmt2 <- rev(hr.stmt)
    if(!stmt.placement %in% c("subtitle", "none")){
        mult <- switch(stmt.placement, bottomleft = 1/50, bottomright = 1 - 1/50, topright = 1 - 1/50)
        algn <- switch(stmt.placement, bottomleft = 4, bottomright = 2, topright = 2)
        incr <- (1:(np - 1))/20
        vert <- switch(stmt.placement, bottomleft = 0 + incr, bottomright = 0 + incr, topright = 1 - incr, none = 0)
        text(x = max(sv$time)*mult, y = vert, labels = hr.stmt2, pos = algn)
    }
    if(stmt.placement=="subtitle"){    
        mtext(hr.stmt2, side = 3, line = 0:(np - 2))
    }
    if(legend.placement!="none") legend(legend.placement, lbls, col = 1:np, lty = 1:np, bty = "n")
}

if(out == "plain"){
    tmp.out <- list(med = medians)
    if(np>1) tmp.out$hr <- hr.txt
} else {
    tmp.out <- list(med = med.stmt)
    if(np>1) tmp.out$hr <- hr.stmt
}

return(tmp.out)

}
