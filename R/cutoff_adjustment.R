



cutoff_adjustment <- function(time, time2, event,
                              outcome, 
                              trim = 2,
                              na.omit = FALSE, 
                              data, 
                              firth = FALSE, reverse = FALSE){
  
  
  
  
  data <- as.data.frame(data)
  
  
  # * NA treatment ----------------------------------------------------------
  
  if(na.omit){
    data <- na.omit(data[ , c(time, time2, event, outcome)])
  }else{
    if(nrow(data) != nrow(na.omit(data[ , c(time, time2, event, outcome)]))){
      stop("missing values are not allowed")
    }
  }
  
  
  
  # * preparation -----------------------------------------------------------
  
  x <- data[, outcome]
  
  
  
  
  ##  Number of events.
  n.events <- sum(data[ , event])
  
  ##  Values of covariate x associated with events.
  meas.event <- sort(x[data[ , event] == 1])
  
  ##  Boundaries for truncation
  bound.l <- head(meas.event, trim)[trim]
  bound.u <- tail(meas.event, trim)[1] 
  
  ##  Truncated set of values of covariates.
  meas.event.trunc <- meas.event[meas.event > bound.l & meas.event < bound.u] 
  
  ##  Set of all possible observations/cutpoints within the
  ##  selection interval.
  cuts <- sort(x[x > bound.l & x < bound.u])
  cuts <- unique(cuts) 
  
  
  
  ##  Number of all cutpoints.
  .n.cuts <- length(cuts)
  
  ##  Obtain p-values, coefficients, s.e. of coefs., and proportion
  ##  epsilon of the sample size <= cutpoint, from Cox regression
  ##  models using all observations/cutpoints within the truncated
  ##  set of values as dichotomous covariates.
  
  .formula <- formula(paste0("Surv(time = ", time, ", time2 = ", time2, ", event = ", event, ") ~ x.tmp"))
  pvals.score <- .threshold <- .ci.lower <- .ci.upper <- pvals <- betas <- .se.betas <- .eps <- rep(NA, times = .n.cuts)
  
  
  
  # * using coxphf (firth penalization) -------------------------------------
  
  
  if(firth){

    if(requireNamespace("coxphf", quietly = TRUE)){
      
      for(i in 1:.n.cuts) {
        .threshold[i] <- cuts[i]
        .x <- factor(x <= cuts[i], labels = c(paste(">", cuts[i], sep = ""), paste("<=", cuts[i], sep = "")))
        
        ##  firth Cox regression model.
        .model <- tryCatch(coxphf::coxphf(formula = Surv(time = , time2 = , event = ), 
                                          data = data.frame(data, x.tmp = .x)), 
                           error = function(e) e, warning = function(w) w)
        
        pvals.score[i] <- NA # score test
        pvals[i] <- ifelse(is(.model, "warning"), NA, summary(.model)$prob) # p.value estiamte
        betas[i] <- ifelse(is(.model, "warning"), NA, coef(.model))
        .ci.lower[i] <- ifelse(is(.model, "warning"), NA, log(summary(.model)$ci.lower))
        .ci.upper[i] <- ifelse(is(.model, "warning"), NA, log(summary(.model)$ci.upper))
        
        .se.betas[i] <- ifelse(is(.model, "warning"), NA, sqrt(vcov(.model)))
        .eps[i] <- ifelse(test = is(.model, "warning"), 
                          yes = NA,
                          no = mean(as.numeric(.x) - 1)) # "-1" due to numeric representation of level 1 = 1, level 2 = 2, ...
      }
      
    }else{
      stop("Package \"coxphf\" needed for this function to work if using Firth penalization. Please install it.",
           call. = FALSE)
    }

    
  }else{
    
    # normal coxph
    
    for (i in 1:.n.cuts) {
      .threshold[i] <- cuts[i]
      .x <- factor(x <= cuts[i], labels = c(paste(">", cuts[i], sep = ""), paste("<=", cuts[i], sep = "")))
      
      ##  Cox regression model.
      .model <- tryCatch(coxph(.formula, data = data.frame(data, x.tmp = .x)), error = function(e) e, warning = function(w) w)
      
      s <- summary(.model)
      
      
      pvals.score[i] <- ifelse(is(.model, "warning"), NA, summary(.model)$sctest["pvalue"]) # score test
      pvals[i] <- ifelse(is(.model, "warning"), NA, summary(.model)$coefficients[5]) # p.value estiamte
      betas[i] <- ifelse(is(.model, "warning"), NA, coef(.model))
      .ci.lower[i] <- ifelse(is(.model, "warning"), NA, (confint(.model)[1]))
      .ci.upper[i] <- ifelse(is(.model, "warning"), NA, (confint(.model)[2]))
      
      .se.betas[i] <- ifelse(is(.model, "warning"), NA, summary(.model)$coefficients[,"se(coef)"])
      .eps[i] <- ifelse(test = is(.model, "warning"),
                        yes = NA,
                        no = mean(as.numeric(.x) - 1)) # "-1" due to numeric representation of level 1 = 1, level 2 = 2, ...
      
    }
  }
  
  
  
  
  # ##  Select cutpoint with lowest p-value. If there are several identical
  # ##  p-values, choose first cutpoint.
  # if (cutpoint == "min") {
  #   .id <- which(pvals == min(na.omit(pvals)))[1]
  # } else if (cutpoint < bound.l) {
  #   stop("Invalid cutpoint specified (below truncated range)!")
  # } else if (cutpoint > bound.u) {
  #   stop("Invalid cutpoint specified (above truncated range)!")
  # } else {
  #   .id <- max(which(cuts <= cutpoint))
  # }
  
  # select minimal p-value
  .id <- which(pvals == min(na.omit(pvals)))[1]
  
  
  
  .p.min <- pvals[.id]
  
  ##  Adjust p-value according to Hilsenbeck & Clark 1996.
  .z <- qnorm(1 - .p.min / 2)
  .l <- length(.eps) - 1
  .D <- rep(NA, times = .l)
  for (j in 1:.l) {
    .s1 <- sqrt(1 - (.eps[j] * (1 - .eps[j+1])) / ((1 - .eps[j]) * .eps[j+1]))
    .s2 <- (.z^2 / 4 - 1) * .s1^3 * 1/6
    .D[j] <- .s1 - .s2
  }
  
  
  p.adj <- .p.min + exp(- .z^2 / 2) / pi * sum(na.omit(.D)) 
  
  ##  Adjust coefficient according to Schumacher et al. 2012, pp. 424f.
  .c.hat <- (betas[.id]^2 - .se.betas[.id]^2) / betas[.id]^2
  beta.adj <- betas[.id] * .c.hat
  
  ##  95% c.i. for adjusted beta, according to Leo.
  ci.beta.adj <- beta.adj + c(-1, 1) * sqrt(.c.hat) * qnorm(p = 0.975) * .se.betas[.id]
  
  
  
  # * reverse outcome -------------------------------------------------------
  
  if(reverse){
    betas <- -betas
    
    .aux <- .ci.lower
    .ci.lower <- -.ci.upper
    .ci.upper <- -.aux
    beta.adj <- -beta.adj
    ci.beta.adj <- -ci.beta.adj
  }
  
  
  
  
  
  # * set return value together ---------------------------------------------
  
  ##  Generate return list.
  ret <- list(n.event = n.events, 
              meas.event = meas.event,
              bound.l = bound.l,
              bound.u = bound.u,
              meas.event.trunc = meas.event.trunc, 
              cutpoints = cuts,
              
              
              beta.unadj = betas,
              beta.unadj.ci.lower = .ci.lower,
              beta.unadj.ci.upper = .ci.upper,
              
              se.beta.unadj = .se.betas,
              
              p.unadj = pvals, 
              p.unadj.min = min(pvals, na.rm = TRUE),
              
              cutoff = cuts[.id],
              p.adj = p.adj, 
              c.hat = .c.hat, # shrinkage factor
              
              beta.adj = beta.adj, 
              ci.beta.adj = ci.beta.adj,
              
              hr.adj = exp(beta.adj),
              ci.hr.adj = sort(exp(ci.beta.adj)),
              
              eps = .eps, 
              pvals.score = pvals.score,
              D = .D)
  
  ret$data.frame <- data.frame(cutpoints = ret$cutpoints,
                               HR = exp(ret$beta.unadj), 
                               CI.lower = exp(ret$beta.unadj.ci.lower),
                               CI.upper = exp(ret$beta.unadj.ci.upper),
                               p.value = ret$p.unadj)
  
  
  return(ret)
  
}
