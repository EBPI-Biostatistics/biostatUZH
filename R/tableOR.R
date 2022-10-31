#' Results table with odds ratios from logistic regression models for binary or
#' ordinal data
#' 
#' Result table with odds ratios, 95\%-CI, test statistics, and p-values from
#' logistic regression models for binary or ordinal variables.
#' 
#' 
#' @param model Object of class \code{"polr"} or \code{"glm"} with family binomial.
#' @param caption Character vector of length 1 containing the table caption.
#' @param label Character vector of length 1 containing the reference label of
#' the Latex table.
#' @param size Character vector of length 1 setting LaTeX font site, e.g. "small",
#' "scriptsize", etc.
#' @param factorNames A character vector of size k number of factors or
#' regressor with custom factor labels.
#' @param table.placement LaTeX table positioning. The default is "ht".
#' See \code{?print.xtable} for more information.
#' @param refLevels A character vector of size k number of regressors with
#' custom reference level names. This argument is usuful for objects of class
#' \code{glm()}, which does not store the reference level labels internally.
#' @param language Language of the table, either "english" (default) or "german".
#' @param short If \code{TRUE}, factor names are removed from factor
#' levels. Default is \code{FALSE}.
#' @param latex If \code{TRUE} (default) LateX output is produced.
#' @param rmStat Logical, if \code{FALSE} (default) output table inludes test
#' statistics.
#' @param wald Logical, if \code{FALSE} (default) Wilson confidence intervals
#' are computed.
#' @return The table as data.frame and, depending on the value of the argument \code{latex},
#' the function a print of the LaTeX code.
#' @author Simon Schwab
#' @seealso \code{\link{tableRegression}}, \code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
#' @examples
#' 
#' data <- carData::TitanicSurvival
#' # relevel: baseline is survived yes.
#' data$survived <- relevel(x = data$survived, ref = "yes") 
#' model <- glm(survived ~ sex + age + passengerClass, data = data, family = binomial())
#' tableOR(model = model, latex = FALSE, short = TRUE,
#'         refLevels = c("female", "1st"), 
#'         caption = "Changes in odds for risk of death in the Titanic tragedy.")
#' 
#' ## using log regression for ordinal data
#' data$passengerClass <- factor(x = data$passengerClass, ordered = TRUE)
#' model <- MASS::polr(passengerClass ~ sex + age, data = data, Hess = TRUE)
#' tableOR(model = model, latex = FALSE, short = TRUE, 
#'         caption = "Changes in odds for being in a lower class, i.e. 2nd or 3rd class")
#'
#' @export
tableOR <- function(model, caption = "", label = "", size = "scriptsize",
                    factorNames = NULL, table.placement = "ht",
                    refLevels = NULL, language = c("english", "german"),
                    short = FALSE, latex = TRUE, rmStat = FALSE, wald = FALSE) {

    if(!(inherits(x = model, what = "polr") ||
         (inherits(x = model, what = "glm") &&  stats::family(model)$family == "binomial")))
        stop("Unknown model. only polr() and glm(..,family = binomial()) models are supported.")
    stopifnot(is.character(caption), length(caption) == 1,
              is.character(label), length(label) == 1,
              is.character(size), length(size) == 1,
              is.null(factorNames) || (is.character(factorNames) && length(factorNames) >= 1),
              is.character(table.placement), length(table.placement) == 1,
              is.null(refLevels) || (is.character(refLevels) && length(refLevels) >= 1),
              !is.null(language))
    language <- match.arg(language)
    stopifnot(is.logical(short), length(short) == 1,
              is.logical(latex), length(latex) == 1,
              is.logical(rmStat), length(rmStat) == 1,
              is.logical(wald), length(wald) == 1)

    
  modelSummary <- summary(model)
  hasCategorial <- !is.null(modelSummary$contrasts)
  
  ## Test which model class we deal with
  isPolr <- class(model)[1] == "polr"
  
  if (nrow(modelSummary$coefficients) < 3) {  # min. requirements: intercept and two factors, or one factor
    # with 2 levels or more (ignoring the baseline level here)
    stop("Table to small, model requires more than one regressor, or a categorial regressor with > 2 levels to create table.")
  }
  
  ## Create table with Odds ratios, 95% CI and test statistic
  if (isPolr) {
    nIntercepts <- length(modelSummary$lev) - 1 # number of intercepts
    k <- length(modelSummary$xlevels) # no. of regressors/factors
    table <- as.data.frame(modelSummary$coefficients[1:(nrow(modelSummary$coefficients) -
                                                        nIntercepts),]) # removing intercepts
    table$`p-value` <- formatPval((1 - stats::pnorm(abs(table$`t value`))) * 2)
    table$OR <- sprintf('%.2f', exp(table$Value))
    if (wald) {
      table$CI <- formatCI(exp(stats::confint.default(model)), text = language)
    } else {
      table$CI <- formatCI(exp(stats::confint(model)), text = language)
    }
    table$`t value` <- sprintf('%.2f', table$`t value`)
  } else { # glm ...
    nIntercepts <- 1
    k <- length(modelSummary$contrasts) # no. of regressors/factors
    table <- as.data.frame(modelSummary$coefficients[-1,]) # removing intercepts
    table$`Pr(>|z|)` <- formatPval(table$`Pr(>|z|)`)
    table$OR <- sprintf('%.2f', exp(table$Estimate))
    if(wald) {
      table$CI <- formatCI(exp(stats::confint.default(model)), text = language)[2:(nrow(table)+1)]
    } else {
      table$CI <- formatCI(exp(stats::confint(model)), text = language)[2:(nrow(table)+1)]
    }
    table$`z value` <- sprintf('%.2f', table$`z value`)
    colnames(table)[4] <- "p-value"
  }
  colnames(table)[6] <- c("95% CI")
  
  ## Label non-categorial variables
  if (hasCategorial) {
    pattern <- names(modelSummary$contrasts)
    table$type <- rep("non-categorial", nrow(table))
    for (i in 1:length(pattern)) {
      table$type[grep(pattern[i], rownames(table))] <- "categorial"
    }
    
    
    ## Improve rownames, remove factor name from level
    if (short) {
        newNames <- rownames(table)
        pattern <- names(modelSummary$contrasts)
        for (i in 1:k) {
            newNames <- sub(pattern[i], "",  newNames)
        }
        rownames(table) <- newNames
    }
    pattern <- names(modelSummary$contrasts)
    levelList <- rownames(modelSummary$coefficients)
    noOfLevels <- rep(NA, k)
    for (i in 1:k) {
        idx <- grep(pattern[i], levelList)
        noOfLevels[i] <- length(idx) + 1
    }

  ## Quick and easy fix if we have not only categorial regressors
  ## we just move these at the bottom of the table.
    table <- rbind(table[table$type == 'categorial',], table[table$type == 'non-categorial',])

  }
  
  ## Get factor names and create new reference level labels myLevelNames
  ## They look as follows: FactorName: ReferenceLevel
    # if factor names are not user specified
    if (is.null(factorNames)) {
        factorNames <- names(modelSummary$contrasts)
    }

    if(hasCategorial){
        if (isPolr) {
            myLevelNames <- rep(NA, k)
            for (i in 1:k) {
                myLevelNames[i] <- paste(factorNames[i], modelSummary$xlevels[[i]][1], sep = ": ")
            } 
        } else {
                                        # unlike polr(), glm() does not contain the reference level labels
                                        # we create dummy ref level names if none are provided
            myLevelNames <- rep(NA, k)
            for (i in 1:k) {
                if (is.null(refLevels)) {
                    myLevelNames[i] = paste(factorNames[i], "reference level", sep = ": ")
                } else {
                    myLevelNames[i] = paste(factorNames[i], refLevels[i], sep = ": ")
                }
            }
        }

    
    ## Add reference levels at the bottom of the table
        new <- data.frame(rep(0, k), rep(NA, k), rep("--", k), rep("--", k), 
                         rep(1, k), rep("--", k), rep('categorial', k),
                         row.names = myLevelNames)
        colnames(new) <- colnames(table)
        table <- rbind(table, new)
        
        
        ## Rearrange table rows so that reference levels are at the top
        ## of each corresponding factor
        startFactor <- nrow(modelSummary$coefficients) - nIntercepts + 1
        
        idx <- NULL
        startLevel <- 1
        for (i in 1:k) {
            idx <- c(idx, startFactor, startLevel:(startLevel + noOfLevels[i] - 2) )
            startFactor <- startFactor + 1
            startLevel <- startLevel + noOfLevels[i] - 1
        }
        table <- rbind(table[idx,], table[-idx,])
        
        ## Determine horizontal lines to group categorial variables
        hline <- c(-1, 0, rep(0, k)) # -1 and 0 is above and below header
        for (i in 1:k) {
            hline[i+2] <- noOfLevels[i] + hline[i+1]
        }
        
                                        # additonal hline at bottom of table
        if (hline[length(hline)] != nrow(table)) {
            hline <- c(hline, nrow(table))
        }
    }
    
    ## Subset table and produce latex code
    if (rmStat) {
        out <- table[,c("OR", "95% CI", "p-value")]
    } else {
        if (isPolr) {
            out <- table[,c("OR", "95% CI", "t value", "p-value")]
        } else {
            out <- table[,c("OR", "95% CI", "z value", "p-value")]
        }
    }
    
    if (latex) {
        xtable::print.xtable(xtable::xtable(out, align = rep("r", ncol(out)+1),
                                            caption = caption, label = label),
                             hline.after = hline, size = size,
              table.placement = table.placement)
        return(invisible(out))
    }
    out
}
