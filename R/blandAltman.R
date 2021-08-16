#' Bland-Altman plot to assess agreement between two continuous variables
#' 
#' Agreement between two continuous measurements is usually assessed by
#' plotting the difference vs. the mean of the two measurements, the so-called
#' Bland-Altman plot. Methods can be used interchangeably if all observations
#' in the Bland-Altman lie between the mean +/- two standard deviations,
#' provided that differences in this range were not clinically important. Note
#' that the function requires pairwise observations.
#' 
#' 
#' @param x Measurements of first method.
#' @param y Measurements of second method. Must have the same length as \code{x}.
#' @param conf.level Confidence level for confidence intervals around the limits
#' of agreement. The default value is 0.95.
#' @param group Factor that enables different colors for points
#' (e.g. to group observations by observation date).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Title of the plot.
#' @param ylim Vector of length 2 containing upper and lower limits of the y-axis.
#' @param plot If \code{TRUE}, the plot is generated.
#' @return A list containing the following elements:
#' \item{difference.mean}{Mean difference.} \item{difference.sd}{Standard
#' deviation of differences.} \item{difference.se}{Standard error of mean
#' differences.} \item{upper.agreement.limit}{Upper limit of agreement.}
#' \item{lower.agreement.limit}{Lower limit of agreement.}
#' \item{agreement.limit.se}{Standard error of limits of agreement.}
#' \item{t.value}{\eqn{t}-quantile used for computation of limits of
#' agreement.}
#' @author Kaspar Rufibach, \email{kaspar.rufibach@@gmail.com}
#' @references Bland, J.M. and Altman, D.G. (1986). Statistical methods for
#' assessing agreement between two methods of clinical measurement.
#' \emph{Lancet}, \bold{1}, 307--310.
#' @keywords dplot aplot htest
#' @aliases BlandAltman blandaltman
#' @examples
#' 
#' n <- 50
#' meas1 <- rnorm(n, 20, 4)
#' meas2 <- meas1 + rnorm(n, 0.5, 1)
#' time <- sample(c(rep(1, n / 2), rep(2, n / 2)))
#' 
#' blandAltman(x = meas1, y = meas2, conf.level = 0.95,
#'             xlab = "mean of measurements", ylab = "difference of measurements", 
#'             main = "Assess agreement between measurement 1 and measurement 2")
#' 
#' blandAltman(x = meas1, y = meas2, conf.level = 0.95, group = time,
#'             xlab = "mean of measurements", ylab = "difference of measurements",
#'             main = "Assess agreement between measurement 1 and measurement 2")
#' 
#' @import stats
#' @export
blandAltman <- function(x, y, conf.level = 0.95, group = NULL,
                        xlab = "average of measurements",
                        ylab = "difference of measurements",
                        main = "", ylim = NULL, plot = TRUE){

    stopifnot(is.numeric(x), length(x) > 0, is.finite(x),
              is.numeric(y), length(y) == length(x), is.finite(y),
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level),  0 < conf.level, conf.level < 1)
    if(!is.null(group)) {
        group <- as.factor(group)
        stopifnot(length(group) == length(x))
    }
    stopifnot(is.character(xlab), length(xlab) == 1,
              is.character(ylab), length(ylab) == 1, 
              is.character(main), length(main) == 1, 
              is.null(ylim) ||
              (is.numeric(ylim) && length(ylim) == 2 && ylim[1] < ylim[2]), 
              is.logical(plot), length(plot) == 1, is.finite(plot))
    
    alpha <- 1 - conf.level
    
    ## use only pairwise complete observations
    ind <- complete.cases(x, y)
    x <- x[ind]
    y <- y[ind]

    difference <- x - y                             # vector of differences
    average <- (x + y) / 2                          # vector of means
    difference.mean <- mean(difference)             # mean difference
    difference.sd <- sd(difference)                 # SD of differences
    al <- qnorm(1 - alpha / 2) * difference.sd
    upper.agreement.limit <- difference.mean + al   # agreement limits
    lower.agreement.limit <- difference.mean - al
    n <- length(difference)                         # number of 'observations'
    
    difference.se <- difference.sd / sqrt(n)        # standard error of the mean
    al.se <- difference.sd * sqrt(3) / sqrt(n)      # standard error of the agreement limit
    tvalue <- qt(1 - alpha / 2, n - 1)              # t value for 95% CI calculation
    difference.mean.ci <- difference.se * tvalue
    al.ci <- al.se * tvalue
    upper.agreement.limit.ci <- c(upper.agreement.limit - al.ci,
                                  upper.agreement.limit + al.ci)
    lower.agreement.limit.ci <- c(lower.agreement.limit - al.ci,
                                  lower.agreement.limit + al.ci)
    
    if (plot) {
        ## The x and the y limits of the plot (xlim, ylim)
        xlim <- range(average)
        if (is.null(ylim)){
            ylim <- range(c(difference, upper.agreement.limit, lower.agreement.limit))
            ylim <- max(abs(ylim)) * c(-1, 1) * 1.1
        }
    
                                        # Plot
        typ <- "p"
        if (!is.null(group)){
            typ <- "n"
            lev <- levels(group)
        }
        plot(x = average, y = difference, type = typ, cex = 1, xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main)
        
                                        # add points if grouping variable is given
        if (!is.null(group)){
            for (i in 1:length(lev))
                points(average[group == lev[i]], difference[group == lev[i]],
                       cex = 0.7, col = i + 1)
            legend("topright", lev, lty = 1, pch = 1, col = 2:(length(lev) + 1), bty = "n")
        }
        abline(h = upper.agreement.limit, lty = "dotted", col = "blue")
        abline(h = difference.mean, col = "blue")
        abline(h = lower.agreement.limit, lty = "dotted", col = "blue")
    }
    
                                        # Return list
    
    out <- list(difference.mean = difference.mean,
                ci.mean = difference.mean + c(-1, 1) * difference.mean.ci, 
                difference.sd = difference.sd, difference.se = difference.se,
                upper.agreement.limit = upper.agreement.limit,
                lower.agreement.limit = lower.agreement.limit, 
                agreement.limit.se = al.se, ci.upper.loa = upper.agreement.limit.ci,
                ci.lower.loa = lower.agreement.limit.ci, t.value = tvalue, n = n)
    if(plot)
        return(invisible(out))
    out
}

#' @export
BlandAltman <- function(x, y, conf.level = 0.95, group = NA,
                        labx = "average of measurements",
                        laby = "difference of measurements",
                        maintit = "", limy = NA, plot = TRUE) {
    .Deprecated("blandAltman")
    blandAltman(x = x, y = y, conf.level = conf.level,
                group = group, xlab = labx, ylab = laby,
                main = maintit, ylim = limy, plot = plot)
}
