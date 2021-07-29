#' Compute arbitrary quantile of a Cumulative Incidence estimate
#' 
#' In survival analysis, often some quantiles (mainly the median) of an
#' estimated cumulative incidence function are of interest. This function
#' computes quantiles of cumulative incidence functions received via
#' \code{cuminc} in the package \pkg{cmprsk}.
#' 
#' 
#' @param time Event times, censored or observed.
#' @param event Event indicator.
#' @param group Quantiles can be computed for several survival curves, defined
#' by \code{group}.
#' @param quant Quantile to be computed. Real number in \eqn{[0, 1]}.
#' @return Dataframe with quantiles estimates for each event-group combination.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @seealso \code{cuminc} and \code{timepoint} in the package \pkg{cmprsk}.
#' @keywords htest survival
#' @examples
#' 
#' if (require("cmprsk")) {
#'   ## illustrate on simulated data
#'   set.seed(1977)
#'   n <- 180
#'   time <- rgamma(n, 2, 1)
#'   status <- factor(sample(rep(0:2, n / 3)), levels = 0:2, labels = 
#'       c("censored", "Time to 2nd tumor", "Death"))
#'   x <- factor(round(runif(n, 1, 2)), levels = 1:2, labels = 
#'       c("Tmt A", "Tmt B"))
#'   
#'   ## plot
#'   plot(cuminc(time, status, x))
#' 
#'   ## median time to event
#'   quantileCumInc(time, status, group = x, quant = 0.25)
#' }
#' 
#' @export
quantileCumInc <- function(time, event, group, quant = 0.5)
{
if (!requireNamespace("cmprsk")) stop("requires cmprsk::cuminc()")

t.low0 <- unlist(lapply(split(time, interaction(group, event, sep = " ")), min))
t.up0 <- unlist(lapply(split(time, interaction(group, event, sep = " ")), max))

est0 <- cmprsk::timepoints(cmprsk::cuminc(time, event, group = group, rho = 0, cencode = 0), times = mean(time))$est
res <- matrix(NA, nrow = nrow(est0), ncol = 1)
rownames(res) <- rownames(est0)
colnames(res) <- paste(quant, "-quantile", sep = "")

for (r in 1:nrow(est0)){
    tmp <- function(x, time, event, group){cmprsk::timepoints(cmprsk::cuminc(time, event, group = group, rho = 0, cencode = 0), times = x)$est[r, 1] - quant}
        
    low <- tmp(x = t.low0[r], time, event, group)
    up <- tmp(x = t.up0[r], time, event, group)

    i2 <- 0
    i4 <- 0
    i5 <- 0

    i1 <- is.na(low) == FALSE
    if (i1 == 1){i2 <- low <= 0}
    i3 <- is.na(up) == FALSE
    if (i3 == 1){i4 <- 0 <= up}
    if (i1 * i3 == 1){i5 <- low < up}
    if (i1 * i2 * i3 * i4 * i5 == 1){res[r, 1] <- uniroot(tmp, interval = c(t.low0[r], t.up0[r]), time, event, group)$root}
}

return(res)
}
