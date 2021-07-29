#' Compute confidence interval for a survival curve at a fixed point
#' 
#' Compute a Wald type confidence interval for a Kaplan-Meier survival curve at
#' a fixed point. The variance is computed according to Peto's formula and the
#' confidence interval is computed using a logit-transformation, to ensure that
#' its bounds lie in \eqn{(0, 1)}. Alternatives are given in the examples
#' below.
#' 
#' 
#' @param time Event times, censored or observed.
#' @param event Censoring indicator, 1 for event, 0 for censored.
#' @param t0 Vector (or single number) of time points to compute confidence
#' interval for.
#' @param conf.level Confidence level for confidence interval.
#' @return \item{t0}{Time points.} \item{S at t0}{Value of survival curve at
#' \code{t0}.} \item{lower.ci}{Lower limits of confidence interval(s).}
#' \item{upper.ci}{Upper limits of confidence interval(s).}
#' @note The variance according to Peto's formula tends to be more conservative
#' than that based on Greenwood's formula.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @keywords htest survival
#' @examples
#' 
#' \dontrun{
#' ## use Acute Myelogenous Leukemia survival data contained in package 'survival'
#' time <- leukemia[, 1]; status <- leukemia[, 2]; x <- as.factor(leukemia[, 3])
#' tmp <- Surv(time, status) ~ 1
#' plot(survfit(tmp, conf.type = "none"), mark = "/", col = 1:2)
#' confIntKM_t0(time, status, t0 = c(10, 25, 50), conf.level = 0.95)
#' 
#' ## an alternative is the log-log confidence interval using Greenwood's
#' ## variance estimate
#' t0 <- 10
#' obj <- survfit(tmp, conf.int = 0.95, conf.type = "log-log", 
#'     type = "kaplan", error = "greenwood")
#' dat <- cbind(obj$time, obj$surv, obj$lower, obj$upper)
#' dat <- dat[dat[, 1] >= t0, ]
#' dat[1, 3:4]
#' 
#' ## this same confidence interval can also be computed using the 
#' ## package km.ci
#' library(km.ci)
#' ci.km <- km.ci(survfit(tmp), conf.level = 0.95, method = "loglog")
#' dat.km <- cbind(ci.km$time, ci.km$surv, ci.km$lower, ci.km$upper)
#' dat.km <- dat.km[dat.km[, 1] >= t0, 3:4]
#' dat.km[1, ]
#' }
#' 
#' @export
confIntKM_t0 <- function(time, event, t0, conf.level = 0.95){

alpha <- 1 - conf.level

# compute ci according to the formulas in Huesler and Zimmermann, Chapter 21
obj <- survfit(Surv(time, event) ~ 1, conf.int = 1 - alpha, conf.type = "plain", type = "kaplan", error = "greenwood", conf.lower = "peto")
St <- summary(obj)$surv
t <- summary(obj)$time
n <- summary(obj)$n.risk
res <- matrix(NA, nrow = length(t0), ncol = 3)

for (i in 1:length(t0)){
    ti <- t0[i]
    if (min(t) > ti){res[i, ] <- c(1, NA, NA)}
    if (min(t) <= ti){    
        if (ti %in% t){res[i, ] <- rep(NA, 3)} else {
            Sti <- min(St[t < ti])
            nti <- min(n[t < ti])        
            Var.peto <- Sti ^ 2 * (1 - Sti) / nti
            Cti <- exp(qnorm(1 - alpha / 2) * sqrt(Var.peto) / (Sti ^ (3 / 2) * (1 - Sti)))
            ci.km <- c(Sti / ((1 - Cti) * Sti + Cti), Cti * Sti / ((Cti - 1) * Sti + 1))
            res[i, ] <- c(Sti, ci.km)}
    }
} # end for

res <- cbind(t0, res)
dimnames(res)[[2]] <- c("t0", "S at t0", "lower.ci", "upper.ci")
return(res)
}
