#' Format coefficient tables of regression models
#'
#' Formats output from the regression model functions: \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' \code{\link[MASS]{glm.nb}}, \code{\link[survival]{coxph}}, and \code{\link{Weibull}}.
#' @param model an object of class \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' \code{negbin} (obtained by \code{\link[MASS]{glm.nb}}), \code{coxph} (obtained by
#' \code{\link[survival]{coxph}}), and list obtained by \code{\link{Weibull}}.
#' @rdname tableRegression
#' @param stats character vector with stats chosen from "estimate", "exp.estimate",
#' "standarderror", "t.value", "ci.95", "p.value".
#' @param col.nam Character vector of same length and order as in \code{stats}.
#' A percentage sign must be escaped by two backslashes.
#' @param row.nam Names of rows, character vector.
#' @param intercept Logical vector of length one indicating whether to provide an intercept
#' or not. If intercept is set TRUE, the first line of the summary output is
#' removed. If the model is a binomial regression, intercept is set FALSE.
#' Intercepts are not available for Weibull or Cox models, because they do not provide any intercept value.
#' @param text Either "english" (default) or "german" indicating the used language names.
#' @param text.ci Either "english", "german" or "none". The language used to denote confidence interval,
#' see \code{\link{displayCI}}.
#' @param eps.pvalue P-values smaller than \code{eps.pvalue} will be formatted as "< eps.pvalue".
#' @param digits Vector of length \code{stats}, digits used for each column.
#' @param big.mark Character vector as in \code{\link[base]{format}}.
#' @param xtable If TRUE a Latex table is returned, otherwise a data.frame is returned. 
#' @param align See \code{\link[xtable]{xtable}}.
#' @param caption See \code{\link[xtable]{xtable}}.
#' @param label See \code{\link[xtable]{xtable}}.
#' @param vars Specify the variables for which regression summaries should
#' be printed.  The argument \code{vars} takes a string vector with
#' the names of the coefficients in the model.
#' @param ... Arguments passed to \code{\link{print.xtable}}.
#' @return Depending on the value of the \code{xtable} argument, the function
#' either prints and returns LaTeX code representing the produced
#' table of coefficients, or it returns the corresponding data frame.
#' @author Sina Rueeger with contributions by Sebastian Meyer.
#' @seealso \code{\link{xtable}}, \code{\link{lm}}, \code{\link{glm}}, \code{\link[MASS]{glm.nb}}
#' \code{\link[survival]{coxph}}, \code{\link{Weibull}}.
#' @details
#' In \code{stats}:
#' \itemize{
#'  \item{}{If \code{t.value} is chosen, the \code{z.value} might be taken, depending on the model.}
#'  \item{}{For lm-models: \code{ci.95} calculates a confidence interval for the estimate.}
#'  \item{}{For glm- and coxph-models: \code{ci.95} calculates a confidence interval for the exp(estimate).}
#' }
#' @examples
#' ## Linear model
#' ## ---------------
#' mod.lm <- lm(Sepal.Length ~ Sepal.Width, data = iris)
#' mod.lm1 <- lm(Sepal.Length ~ .^2, data = iris) 
#' 
#' tableRegression(model = mod.lm)
#' 
#' ## choosing columns, columns and row naming in german
#' tableRegression(model = mod.lm1, stats = c("estimate", "t.value", "p.value"),
#'                 text = "german")
#' 
#' ## adapt row names, plus special format for ci
#' tableRegression(model = mod.lm, row.nam = c("Intercept", "Width Sepal"),
#'                 text.ci = "none")
#' 
#' ## Poisson model
#' ## (example from ?glm)
#' ## --------------
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' d.AD <- data.frame(treatment, outcome, counts)
#' mod.glm.pois <- glm(counts ~ outcome + treatment, family=poisson())
#' tableRegression(model = mod.glm.pois)
#' 
#' 
#' ## Negative binomial model
#' ## --------------
#' if (require("MASS")) {
#'     mod.glm.nb <- glm.nb(Days ~ Sex + Age, data = quine)
#'     tableRegression(
#'         mod.glm.nb,
#'         caption = paste("NegBin model. Estimated dispersion:",
#'             sprintf("%4.2f ($se=%4.2f$).", mod.glm.nb$theta, mod.glm.nb$SE.theta)),
#'         label = "tab:glm.nb"
#'     )
#' }
#' 
#' 
#' ## Logistic model
#' ## -------------
#' dat <- survival::rats
#' dat$rx <- factor(dat$rx, labels = c(" (A)", " (B)"))
#' mod.glm.bin <- glm(status ~ litter + rx, family = binomial, data = dat)
#' 
#' tableRegression(model = mod.glm.bin,
#'                 stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
#'                 text = "english", digits = rep(3, 5),
#'                 caption = "Here goes the caption.", label = "mod:logit")
#' 
#' ## including intercept
#' tableRegression(model = mod.glm.bin,
#'                 stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
#'                 text = "english", digits = rep(3, 5),
#'                 caption = "Here goes the caption.", label = "mod:logit",
#'                 intercept = TRUE)
#' 
#' 
#' ## Cox model
#' ## (example from ?survival::coxph)
#' ## -------------
#' dat <- list(time=c(4,3,1,1,2,2,3), 
#'             status=c(1,1,1,0,1,1,0), 
#'             x=c(0,2,1,1,1,0,0), 
#'             sex=c(0,0,0,0,1,1,1)) 
#' 
#' library("survival")
#' mod.cox <- coxph(Surv(time, status) ~ x, dat)
#' mod.cox1 <- coxph(Surv(time, status) ~ x + factor(sex), dat)
#' mod.cox2 <- coxph(Surv(time, status) ~ x + strata(sex), dat)
#' 
#' tableRegression(model = mod.cox)
#' tableRegression(model = mod.cox1)
#' tableRegression(model = mod.cox2)
#' 
#' 
#' ## Weibull
#' ## (example from biostatUZH::WeibullReg)
#' ## -------------
#' data("larynx")
#' mod.wb <- WeibullReg(Surv(time, death) ~ factor(stage) + age, data=larynx)
#' tableRegression(model = mod.wb)
#' @export
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
                            vars = NULL,
                            ...
)
{
    
    raw.col.nam.german <- c("Koeffizient", "Exp(Koeffizient)", "Standardfehler", "$t$-Wert", "95\\%-Konfidenzintervall", "$p$-Wert")
    
    raw.col.nam.english <- c("Coefficient", "Exp(Coefficient)", "Standarderror", "$t$-value", "95\\%-confidence interval", "$p$-value")
    
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
    
    # linear model -----------------------------------------------------------
    ## linear model
    
    if (clm == "lm")
    {
        estimate <- summary(model)$coef[,1]
        exp.estimate <- exp(estimate)
        standarderror <- summary(model)$coef[,2]
        t.value <- summary(model)$coef[,3]
        p.value <- summary(model)$coef[,4]
        ci.95 <- formatCI(confint(model), digits = digits.ci, text = text.ci)
    }
    
    
    
    
    # negbin ------------------------------------------------------------------
    ## negbin
    
    
    if (clm %in% c("negbin")) {
        
        estimate <- summary(model)$coef[, 1]
        exp.estimate <- exp(estimate)
        standarderror <- summary(model)$coef[, 2]
        t.value <- summary(model)$coef[, 3]
        p.value <- summary(model)$coef[, 4]
        
        # change colnames
        col.nam[stats == "exp.estimate"] <- "Rate Ratio"
        
        ## confint for exp.estimate (actually depends on MASS:::confint.glm)
        ci.95 <- if (requireNamespace("MASS", quietly = FALSE)) {
            formatCI(exp(confint(model)), digits = digits.ci, text = text.ci)
        } else {
            rep.int(NA_character_, length(estimate))
        }
        
    }
    
    
    # glm / gee ---------------------------------------------------------------
    ## glm / gee
    
    
    if (clm %in% c("glm", "geeglm")) { # both using the same family
        
        back.trafo <- switch(paste0(c(model$family$family, model$family$link), collapse = "."),
                             "binomial.logit" = list("fct" = exp, "name" = "Odds Ratio"),
                             "quasibinomial.logit" = list("fct" = exp, "name" = "Odds Ratio"),
                             
                             "gaussian.identity" = list("fct" = identity, "name" = "Estimate"),
                             "quasi.identity" = list("fct" = identity, "name" = "Estimate"),
                             
                             # "Gamma.inverse" = list("fct" = function(x) 1/x, "name" = "inv(Coefficient)"),
                             # "inverse.gaussian.1/mu^2" = list("fct" = function(x) sqrt(1/x), "name" = "sqrt(inv(Coefficient))"),
                             
                             "poisson.log" = list("fct" = exp, "name" = "Rate Ratio"),
                             "quasipoisson.log" = list("fct" = exp, "name" = "Rate Ratio"))
        
        
        
        
        # change colnames
        col.nam[stats == "exp.estimate"] <- c(back.trafo$name)
        
        # calculating
        estimate <- summary(model)$coef[, 1]
        exp.estimate <- back.trafo$fct(estimate) # bad name for variable
        standarderror <- summary(model)$coef[, 2]
        t.value <- summary(model)$coef[, 3]
        p.value <- summary(model)$coef[, 4]
        
        
        # CI calculation
        ci.95 <- if (requireNamespace("MASS", quietly = FALSE)) {
            if (clm %in% c("glm")) {
                formatCI(back.trafo$fct(confint(model)), digits = digits.ci, 
                         text = text.ci)
            }else {
                formatCI(back.trafo$fct(confint.geeglm.broom(model)), digits = digits.ci, 
                         text = text.ci)
            }
        }else {
            rep.int(NA_character_, length(estimate))
        }
        
    }
    
    
    # coxmod ------------------------------------------------------------------
    ## coxmod
    
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
    
    
    # weibull -----------------------------------------------------------------
    ## weibull
    
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
    
    
    # bring everything together -----------------------------------------------
    ## bring everything together
    
    output <- data.frame(estimate, exp.estimate, standarderror, t.value, ci.95, p.value,
                         stringsAsFactors = FALSE)
    
    if(!intercept & !(clm %in% c("list", "coxph"))) ## in weibull and coxph there is anyway no intercept plotted #
    {
        output <- output[-1,]
    }
    
    
    # extrahieren des return outputs ------------------------------------------
    ## extrahieren des return outputs
    
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
    
    
    # formatieren des outputs -------------------------------------------------
    ## formatieren des outputs
    
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
    
    
    # table -------------------------------------------------------------------
    ## table
    
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

