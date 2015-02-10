WeibullDiag <-
function (formula, data = parent.frame(), labels = rownames(np$X)) 
{
    if (!requireNamespace("prodlim")) stop("requires prodlim::prodlim()")
    np <- prodlim::prodlim(formula, data)
    cols <- NA
    for (i in 1:length(np$size.strata)) {
        cols <- c(cols, rep(i, np$size.strata[i]))
    }
    cols <- cols[!is.na(cols)]
    window <- c(-0.5, 0.5)
    y.range <- range(log(-log(np$surv)), finite = TRUE)
    x.range <- log(range(np$time))
    plot(0, 0, type = "n", ylim = y.range + window, xlim = x.range + 
        window, xlab = "Log Survival Time", ylab = "Log Cumulative Hazard", 
        main = "Weibull Regression\nDiagnostic Plot")
    for (i in 1:length(np$size.strata)) {
        t <- np$time[cols == i]
        s <- np$surv[cols == i]
        h <- np$haz[cols == i]
        points(log(t), log(-log(s)), col = i, pch = i, lty = i, 
            lwd = 2, type = "b")
    }
    legend("topleft", legend = labels, col = 1:length(np$size.strata), 
        pch = 1:length(np$size.strata), lwd = 2, lty = i:length(np$size.strata))
    pCH <- list(x = log(np$time), y = log(-log(np$surv)), strata = cols, 
        stata.def = np$X)
}
