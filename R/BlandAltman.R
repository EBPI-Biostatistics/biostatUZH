BlandAltman <- function(x, y, conf.level = 0.95, group = NA, labx = "average of measurements", 
laby = "difference of measurements", maintit = "", limy = NA, plot = TRUE){

    alpha <- 1 - conf.level
    
    ## use only pairwise complete observations
    ind <- complete.cases(x, y)
    x <- x[ind]
    y <- y[ind]

    difference <- x - y                                 # vector of differences
    average <- (x + y) / 2                              # vector of means
    difference.mean <- mean(difference)                 # mean difference
    difference.sd <- sd(difference)                     # SD of differences
    al <- qnorm(1 - alpha / 2) * difference.sd
    upper.agreement.limit <- difference.mean + al       # agreement limits
    lower.agreement.limit <- difference.mean - al
    n <- length(difference)                             # number of 'observations'
    
    difference.se <- difference.sd / sqrt(n)            # standard error of the mean
    al.se <- difference.sd * sqrt(3) / sqrt(n)          # standard error of the agreement limit
    tvalue <- qt(1 - alpha / 2, n - 1)                  # t value for 95% CI calculation
    difference.mean.ci <- difference.se * tvalue
    al.ci <- al.se * tvalue
    upper.agreement.limit.ci <- c(upper.agreement.limit - al.ci, upper.agreement.limit + al.ci)
    lower.agreement.limit.ci <- c(lower.agreement.limit - al.ci, lower.agreement.limit + al.ci)
    
    if (identical(plot, TRUE)){
    # The x and the y limits of the plot (limx, limy)
    limx <- range(average)
    if (max(is.na(limy)) == 1){
        limy <- range(c(difference, upper.agreement.limit, lower.agreement.limit))
        limy <- max(abs(limy)) * c(-1, 1) * 1.1
        }
    
    # Plot
    typ <- "p"
    if (is.na(group[1]) == FALSE){
        typ <- "n"
        lev <- levels(as.factor(group))
    }
    plot(average, difference, type = typ, cex = 1, xlim = limx, ylim = limy, xlab = labx, ylab = laby, main = maintit)

    # add points if grouping variable is given
    if (is.na(group[1]) == FALSE){
        for (i in 1:length(lev)){points(average[group == lev[i]], difference[group == lev[i]], cex = 0.7, col = i + 1)}
        legend("topright", lev, lty = 1, pch = 1, col = 2:(length(lev) + 1), bty = "n")    
        }
    abline(h = upper.agreement.limit, lty = "dotted", col = "blue")
    abline(h = difference.mean, col = "blue")
    abline(h = lower.agreement.limit, lty = "dotted", col = "blue")
    }
    
    # Return list
    ba <- list(difference.mean = difference.mean, ci.mean = difference.mean + c(-1, 1) * difference.mean.ci, 
    difference.sd = difference.sd, difference.se = difference.se,
    upper.agreement.limit = upper.agreement.limit, lower.agreement.limit = lower.agreement.limit, 
    agreement.limit.se = al.se, ci.upper.loa = upper.agreement.limit.ci, ci.lower.loa = lower.agreement.limit.ci, 
    t.value = tvalue, n = n)
    return(ba)
}

