computeICCrater <- function(dat)
{
    ## we compute the ICC(2, 1): Each subject is measured by each rater, 
    ## and raters are considered representative of a larger population of 
    ## similar raters. Reliability calculated from a single measurement.
    if (!requireNamespace("lme4")) stop("requires lme4::lmer()")
    
    ## compute ICC using a linear mixed model approach
    fm.con <- lme4::lmer(score ~ 1 + (1|pat) + (1|rater), data = dat, REML = TRUE)
    v <- lme4::VarCorr(fm.con)
    sig.res <- as.numeric(attr(v, "sc")) ^ 2
    v <- unlist(lapply(v, as.numeric))
    sig.pat <- v[c("pat")]
    sig.rat <- v[c("rater")]
    icc21 <- as.numeric(sig.pat / (sig.pat + sig.rat + sig.res))    
    icc31 <- as.numeric(sig.pat / (sig.pat + sig.res))
    res.mle <- c("sigmas" = c(sig.pat, sig.rat, sig.res, "icc(2, 1)" = icc21, "icc(3, 1)" = icc31))

    ## generate output    
    res <- rbind(res.mle)
    dimnames(res)[[2]] <- c("sig.pat", "sig.rater", "sig.res", "icc(2, 1)", "icc(3, 1)")

    return(res)
}

