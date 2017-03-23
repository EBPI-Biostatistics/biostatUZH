################################################################################
### Part of the R package "biostatUZH".
### Free software under the terms of the GNU General Public License (version 2
### or later) a copy of which is available at http://www.R-project.org/Licenses
###
### Copyright (C) 2012-2013 Sina Ruegger, 2015 Sebastian Meyer, 2017 Leonhard Held
### Time-stamp: <[tableRegression.R] 2015-05-08 11:07 (CEST) by SM>
################################################################################


tableRegression <- function(model,
                            stats = NULL,
                            col.nam = NULL,
                            row.nam = NULL,
                            intercept = NULL,
                            text = "english", 
                            text.ci = text, 
                            eps.pvalue = 0.0001,
                            digits = NULL,
                            big.mark = "'",
                            xtable = TRUE,
                            align = NULL,
                            caption = NULL,
                            label = NULL,
                            ...
                            )
{
    
    raw.col.nam.german <- c("Koeffizient", "Exp(Koeffizient)", "Standardfehler", "$t$-Wert", "95\\%-Konfidenzintervall", "$p$-Wert")

    raw.col.nam.english <- c("Coefficent", "Exp(Coefficient)", "Standarderror", "$t$-value", "95\\%-confidence interval", "$p$-value")
    
    raw.stats <- c("estimate", "exp.estimate", "standarderror", "t.value", "ci.95", "p.value")
    
    clm <- class(model)[1]
    #if(clm == "glm")
    if (clm %in% c("glm", "geeglm"))
    {
        cl <- model$family$family
    }
    else
    {
        cl <- clm
    }
    ## lm >> linear model
    ## binomial >> generalized linear model
    ## poisson >> generalized linear model
    ## list >> weibull
    ## coxph >> survival
    ## negbin >> negative binomial model (fit by MASS::glm.nb)
    
    ## LM
    ## -------------
    if (clm == "lm")
    {
        ## intercept
        if(is.null(intercept)) intercept <- TRUE

        ## stats
        k <- c(1, 5, 6)
        if(is.null(stats)) stats <- raw.stats[k]

        ## col.nam >> dependent on stats & text
        ind <- sapply(stats, function(x) which(x == raw.stats)) #which(raw.stats %in% stats)
        if(is.null(col.nam))
            {
                if(text == "german") col.nam <- raw.col.nam.german[ind]
                if(text == "english") col.nam <- raw.col.nam.english[ind]
            }
        
        ## row.nam >> dependent on intercept & text
        if(is.null(row.nam))
            {
                row.nam <- names(model$coef)[-1]
                if(intercept)
                {
                    if(text == "german") intercept.nam <- "Achsenabschnitt"
                    if(text == "english") intercept.nam <- "Intercept"
                    row.nam <- c(intercept.nam, row.nam)
                }
                
            }
        
        ## digits
        if(is.null(digits)) digits <- rep(2, length(stats))

    } else {
    ## rest
    ## -------------
        
        ## intercept
        if(is.null(intercept))
        {
            ## intercept is omitted when having weibull
            ## or logit or Poisson regression or coxph
            intercept <- ! ((clm %in% c("list", "coxph")) | (cl %in% c("binomial", "quasibinomial", "poisson", "quasipoisson")))
        }
        
        ## stats
        k <- c(2, 5, 6)
        if(is.null(stats)) stats <- raw.stats[k]

        ## col.nam >> dependent on stats & text
        ind <- sapply(stats, function(x) which(x == raw.stats))#which(raw.stats %in% stats)
        if(is.null(col.nam))
        {
            exp.nam <- switch(
                cl,
                "poisson" =, "quasipoisson" =, "negbin" = "Rate Ratio",
                "binomial" =, "quasibinomial" = "Odds Ratio",
                "coxph" = "Hazard Ratio",
                "Exp(Coefficient)"
            )
            
            ind.exp <- raw.stats == "exp.estimate"
            raw.col.nam.german[ind.exp] <- exp.nam
            raw.col.nam.english[ind.exp] <- exp.nam
            
            if(text == "german") col.nam <- raw.col.nam.german[ind]
            if(text == "english") col.nam <- raw.col.nam.english[ind]
        }
        
        ## row.nam >> dependent on intercept & text
        if(is.null(row.nam))
        {
            if(clm == "list")
            {
                row.nam <- rownames(model$coef)[-c(1,2)]
            }else{
                if(clm == "coxph")
                {
                    row.nam <- names(model$coefficients)
                }else{
                    row.nam <- names(model$coef)[-1]
                }
            }
            
            if(intercept)
            {
                if(text == "german") intercept.nam <- "Achsenabschnitt"
                if(text == "english") intercept.nam <- "Intercept"
                row.nam <- c(intercept.nam, row.nam)
            }
        }
        
        ## digits
        if(is.null(digits)) digits <- rep(2, length(stats))

    }


    ## warning intercept
    if(intercept & clm %in% c("list", "coxph"))
    {
        warning("Weibull and Cox models do not include an intercept. Set intercept = FALSE")
    }
    
    ## digitis ci
    if("ci.95" %in% stats)
    {
        digits.ci <- digits[stats %in% "ci.95"]
    }else{
        digits.ci <- 2
    }
            
    ## text.ci
    ## if(is.null(text.ci)) text.ci <- text
    
    #col.nam <- sub("%", "\\\\%", col.nam)

       
    ## linear model
    ## ----------------------------
    if (clm == "lm")
    {
        estimate <- summary(model)$coef[,1]
        exp.estimate <- exp(estimate)
        standarderror <- summary(model)$coef[,2]
        t.value <- summary(model)$coef[,3]
        p.value <- summary(model)$coef[,4]
        ci.95 <- formatCI(confint(model), digits = digits.ci, text = text.ci)
    }

    ## glm
    ## ----------------------------
    if ((cl %in% c("binomial", "poisson", "quasibinomial", "quasipoisson")) | (clm == "negbin"))
    {
        estimate <- summary(model)$coef[,1]
        exp.estimate <- exp(estimate)
        standarderror <- summary(model)$coef[,2]
        t.value <- summary(model)$coef[,3]
        p.value <- summary(model)$coef[,4]
        if (clm %in% c("glm", "negbin")){
        ## confint for exp.estimate (actually depends on MASS:::confint.glm)
            ci.95 <- if (requireNamespace("MASS", quietly = FALSE)) {
                         formatCI(exp(confint(model)), digits = digits.ci, text = text.ci)
                     } else {
                         rep.int(NA_character_, length(estimate))
                     }
        }
        if (clm %in% c("geeglm")){
            lower <- exp(estimate - 1.96*standarderror)
            upper <- exp(estimate + 1.96*standarderror)
            ci.95 <- formatCI(cbind(lower, upper), digits = digits.ci, text = text.ci)
        }
        
    }


    ## coxmod
    ## ----------------------------
    if (clm == "coxph")
    {
        estimate <- summary(model)$coefficients[,1]
        exp.estimate <- exp(estimate)
        standarderror <- summary(model)$coefficients[,3]
        t.value <- summary(model)$coefficients[,4]
        p.value <- summary(model)$coefficients[,5]
        ci.95 <- formatCI(cbind(summary(model)$conf.int[,3], summary(model)$conf.int[,4]), digits = digits.ci, text = text.ci) 

       # col.nam[stats == "exp.estimate"] <- c("Hazard Ratio")
       ## cl.2 <- "survival"
    }

    ## weibull
    ## ----------------------------
    if (clm == "list")
    {
        estimate <- model$coef[-c(1:2),1]
        exp.estimate <- exp(estimate)
        standarderror <- model$coef[-c(1:2), 2]
        t.value <- NA
        p.value <- 2 * pnorm(- abs(estimate / standarderror))
        ci1 <- estimate - qnorm(0.975) * standarderror
        ci2 <- estimate + qnorm(0.975) * standarderror
        ci.95 <- formatCI(exp(cbind(ci1, ci2)), digits = digits.ci, text = text.ci) 

        col.nam[stats == "exp.estimate"] <- c("Hazard Ratio")
       ## cl.2 <- "survival"
    }

    ## bring everything together
    ## --------------------------
    output <- data.frame(estimate, exp.estimate, standarderror, t.value, ci.95, p.value,
                         stringsAsFactors = FALSE)

    if(!intercept & !(clm %in% c("list", "coxph"))) ## in weibull and coxph there is anyway no intercept plotted #
    {
        output <- output[-1,]
    }
    
    ## extrahieren des return outputs
    ## --------------------------------
    if (nrow(output) > 1)
    {
        output.return <- output[, ind]
        colnames(output.return) <- col.nam
        rownames(output.return) <- row.nam
    }else{
        output.return <- data.frame(output[,ind])
        names(output.return) <- col.nam
        rownames(output.return) <- row.nam
    }

    ## formatieren des outputs
    ## ------------------------
    for (i in 1:ncol(output.return))
    {
        if(stats[i] == "p.value")
        {
            output.return[,i] <- biostatUZH::formatPval(as.numeric(as.character(output.return[,i])), break.eps = eps.pvalue)# dig[i])
        }else{
            if(stats[i] != "ci.95")
            {
                output.return[,i] <-  sapply(output.return[,i], function(x)
                                             format(as.numeric(as.character(x)), big.mark = big.mark, digits = digits[i], nsmall =
                                            digits[i], scientific = FALSE))

               
                
                if(stats[i] == "exp.estimate" & nrow(output.return) > 1)
                {
                    output.return[-1,i] <- sapply(output.return[-1,i], function(x) format(as.numeric(as.character(x)), digits = digits[i], big.mark = big.mark, nsmall =
                                            digits[i], scientific = FALSE))
                } # end if
            } # end if
        } # end ifelse
    } # end for

    #if ("exp.estimate" %in% stats  & cl %in% c("weibull", "coxph") & intercept)
    #{
    #    if ("exp.estimate" %in% stats)
    #    is.na(output.return[1,"exp.estimate" == stats]) <- TRUE#

    #    if ("ci.95" %in% stats)
    #    is.na(output.return[1,"ci.95" == stats]) <- TRUE
    #}


    ## table
    ## --------
    if (xtable && requireNamespace("xtable")) {
        if (is.null(align))
            align <- paste(rep("r", length(stats)+1), collapse = "")
        xtab <- xtable::xtable(output.return, caption = caption, label = label, align = align)
        ## options for print.xtable (can be overridden by ... arguments)
        oopt <- options(
            xtable.include.rownames = TRUE,
            xtable.floating = TRUE,
            xtable.type = "latex",
            xtable.size = "footnotesize",
            xtable.table.placement = "!h",
            xtable.sanitize.colnames.function = identity  # do not escape "$"
        )
        on.exit(options(oopt))
        print(xtab, ...)
    } else {
        output.return
    }
}

