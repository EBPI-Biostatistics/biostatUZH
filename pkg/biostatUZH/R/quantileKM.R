quantileKM <- function(time, event, group = NA, quant = 0.5, conf.level = 0.95, conftype = c("log","log-log","plain","none")[2],
conflower = c("usual", "peto", "modified")[1]){

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
