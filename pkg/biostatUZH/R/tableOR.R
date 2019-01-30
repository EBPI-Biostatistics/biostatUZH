## Part of the R package "biostatUZH".
## by Simon Schwab, 2018
## 
tableOR = function(model, caption="", label="", size="scriptsize", factorNames=NULL,
                   table.placement = "ht", refLevels=NULL, lang="english", short = FALSE,
                   latex=TRUE, rmStat=FALSE) {
  
  mySummary = summary(model)
  hasCategorial = !is.null(mySummary$contrasts)
  
  ## Test which model class we deal with
  isPolr = class(model)[1] == "polr"
  isGlm = class(model)[1] == "glm"
  
  error = "Unknown model, fit must be from polr() or glm(..,family = binomial())."
  if (!isPolr & !isGlm) {
    stop(error)
  } else if (isGlm) {
    if (mySummary$family$family != "binomial") {stop(error)}
  } 
  
  # if (length(attr(mySummary$terms, 'dataClasses')) <= 2) {
  if (nrow(mySummary$coefficients) < 3) {  # min. requirements: intercept and two factors, or one factor
    # with 2 levels or more (ignoring the baseline level here)
    stop("Table to small, model requires more than one regressor, or a categorial regressor with >2
         levels to create table.")
  }
  
  ## Create table with Odds ratios, 95% CI and test statistic
  if (isPolr) {
    nIntercepts = length(mySummary$lev) - 1 # number of intercepts
    k = length(mySummary$xlevels) # no. of regressors/factors
    table = as.data.frame(mySummary$coefficients[1:(nrow(mySummary$coefficients)-nIntercepts),]) # removing intercepts
    table$`p-value` = biostatUZH::formatPval((1 - pnorm(abs(table$`t value`))) * 2) # Check with Leo
    table$OR  = sprintf('%.2f', exp(table$Value))
    table$CI = formatCI(exp(confint.default(model)), text = lang)
    #table$CI = formatCI(exp(confint(model)), text = lang)
    table$`t value` = sprintf('%.2f', table$`t value`)
    
  } else if (isGlm) {
    nIntercepts = 1
    k = length(mySummary$contrasts) # no. of regressors/factors
    table = as.data.frame(mySummary$coefficients[-1,]) # removing intercepts
    table$`Pr(>|z|)` = biostatUZH::formatPval(table$`Pr(>|z|)`)
    table$OR  = sprintf('%.2f', exp(table$Estimate))
    table$CI = formatCI(exp(confint.default(model)), text = lang)[2:(nrow(table)+1)]
    #table$CI = formatCI(exp(confint(model)), text = lang)[2:(nrow(table)+1)]
    table$`z value` = sprintf('%.2f', table$`z value`)
    colnames(table)[4] = c("p-value")
  }
  colnames(table)[6] = c("95% CI")
  
  ## Label non-categorial variables
  if (hasCategorial) {
    pattern = names(mySummary$contrasts)
    table$type = rep('non-categorial', nrow(table))
    for (i in 1:length(pattern)) {
      table$type[grep(pattern[i], rownames(table))] = 'categorial'
    }
  }
  
  ## Improve rownames, remove factor name from level
  if (hasCategorial & short) {
    newNames = rownames(table)
    pattern = names(mySummary$contrasts)
    for (i in 1:k) {
      newNames = sub(pattern[i], "",  newNames)
    }
    rownames(table) = newNames
  }
  
  ## Determine no. of levels per factor
  if (hasCategorial) {
    pattern = names(mySummary$contrasts)
    levelList = rownames(mySummary$coefficients)
    noOfLevels = rep(NA, k)
    for (i in 1:k) {
      idx = grep(pattern[i], levelList)
      noOfLevels[i] = length(idx) + 1
    }
  }
  
  ## Quick and easy fix if we have not only categorial regressors
  ## we just move these at the bottom of the table.
  if (hasCategorial) {
    table = rbind(table[table$type == 'categorial',], table[table$type == 'non-categorial',])
  }
  ## Get factor names and create new reference level labels myLevelNames
  ## They look as follows: FactorName: ReferenceLevel
  
  # if factor names are not user specified
  if (is.null(factorNames)) {factorNames = names(mySummary$contrasts)}
  
  if (isPolr & hasCategorial) {
    myLevelNames = rep(NA, k)
    for (i in 1:k) {
      myLevelNames[i] = paste(factorNames[i], mySummary$xlevels[[i]][1], sep = ": ")
    } 
    
  } else if (isGlm & hasCategorial) {
    # unlike polr(), glm() does not contain the reference level labels
    # we create dummy ref level names if none are provided
    myLevelNames = rep(NA, k)
    for (i in 1:k) {
      if (is.null(refLevels)) {
        myLevelNames[i] = paste(factorNames[i], "reference level", sep = ": ")
      } else {
        myLevelNames[i] = paste(factorNames[i], refLevels[i], sep = ": ")
      }
    }
  }
  
  ## Add reference levels at the bottom of the table
  if (hasCategorial) {
    new = data.frame(rep(0,k), rep(NA,k), rep("--",k), rep("--",k), 
                     rep(1,k), rep("--",k), rep('categorial',k),
                     row.names = myLevelNames)
    colnames(new) = colnames(table)
    table = rbind(table, new)
    
    
    ## Rearrange table rows so that reference levels are at the top
    ## of each corresponding factor
    startFactor = nrow(mySummary$coefficients) - nIntercepts + 1
    
    idx = NULL
    startLevel = 1
    for (i in 1:k) {
      idx = c(idx, startFactor, startLevel:(startLevel + noOfLevels[i] - 2) )
      startFactor = startFactor + 1
      startLevel = startLevel + noOfLevels[i] - 1
    }
    table = rbind(table[idx,], table[-idx,])
    
    ## Determine horizontal lines to group categorial variables
    hline = c(-1, 0, rep(0, k)) # -1 and 0 is above and below header
    for (i in 1:k) {
      hline[i+2] = noOfLevels[i] + hline[i+1]
    }
    
    # additonal hline at bottom of table
    if (hline[length(hline)] != nrow(table)) {
      hline = c(hline, nrow(table))
    }
  }
  
  ## Subset table and produce latex code
  if (rmStat) {
    out = table[,c("OR", "95% CI", "p-value")]
  } else {
    if (isPolr) {out = table[,c("OR", "95% CI", "t value", "p-value")]}
    else if (isGlm) {out = table[,c("OR", "95% CI", "z value", "p-value")]}
  }
  
  if (latex) {
    print(xtable::xtable(out, align = rep("r", ncol(out)+1), caption = caption, label = label),
          hline.after = hline, size = size, table.placement = table.placement)
    
  } else {return(out)}
}