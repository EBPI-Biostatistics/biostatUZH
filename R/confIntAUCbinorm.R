confIntAUCbinorm <- function(cases, controls, conf.level = 0.95, replicates =
                             1000, grid = 100, ci.method = c("boot", "wald"),
                             var.equal = TRUE){ 

    ## decide CI method
    ci.method <- match.arg(ci.method)
    
    alpha <- 1 - conf.level

    ## CI with bootstrap
    ## -------------------
    if (ci.method == "boot")
    {
        ## compute bootstrap CI for binormal AUC
        ## data in two columns: y (measurements), status (0 = control, 1 = case)
        
        ts <- seq(0, 1, by = 1 / grid)		
        a <- (mean(cases) - mean(controls)) / sd(cases)
        b <- sd(controls) / sd(cases)	
        biroc <- pnorm(a + b * qnorm(ts))

        dat <- cbind(c(controls, cases), c(rep(0, length(controls)), rep(1, length(cases))))

        ## compute bootstrap ci
        res <- boot(dat, AUCbinorm, R = replicates)
        aucbin <- boot.ci(res,type = "bca")

        return(list("a" = a, "b" = b, "x.val" = ts, "y.val" = biroc,
                    "auc" = res$t0, "lowBootCI" = aucbin$bca[4],
                    "upBootCI" = aucbin$bca[5]))  
    }

    ## CI 
    ## -------------------
    if(ci.method == "wald")
    {
        n0 <- length(controls)
        n1 <- length(cases)
        
        mu0 <- mean(controls)
        mu1 <- mean(cases)
        
        if(var.equal)
        {
            s0 <- s1 <- sd(c(cases, controls))
            a.se <- sqrt(1/n0 + 1/n1)
            auc.se <- a.se
        } else {
            stop("var.equal = FALSE not implemented for ci.method = 'wald'")
            
            ## s0 <- sd(controls)
            ## s1 <- sd(cases)
        }

        a <- (mu1 - mu0)/s1
        b <- s0/s1
        auc <- pnorm(a/(sqrt(1 + b^2)))
        lowCI <- pnorm((a - qnorm(1-alpha/2) * auc.se) / (sqrt(1 + b^2)))
        upCI <- pnorm((a + qnorm(1-alpha/2) * auc.se) / (sqrt(1 + b^2)))
        
        return(list("a" = a, "b" = b, "auc" = auc, "lowCI" = lowCI, "upCI" = upCI))
    }


}
