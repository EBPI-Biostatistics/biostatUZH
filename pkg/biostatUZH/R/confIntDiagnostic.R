confIntDiagnostic <- function(tp, fp, tn, fn, conf.level = 0.95, cohort=FALSE, pr=NA)
{
    stopifnot(is.wholenumber(tp), is.wholenumber(fp),
              is.wholenumber(tn), is.wholenumber(fn),  conf.level<1,
              conf.level>0)
    res <- data.frame(matrix(NA, nrow = 7, ncol = 4))

    colnames(res) <- c("type", "lower", "estimate", "upper")

    res[1, 2:4] <- wilson(x=tp, n=tp+fn, conf.level = conf.level)
    res[2, 2:4] <- wilson(x=tn, n=tn+fp, conf.level = conf.level)
    LRplus <- confIntRiskRatio(x=c(tp,fp), n=c(tp+fn, fp+tn), conf.level = conf.level)
    LRminus <- confIntRiskRatio(x=c(fn,tn), n=c(tp+fn, tn+fp), conf.level = conf.level)
    res[3, 2:4] <- LRplus
    res[4, 2:4] <- LRminus
    
    if(!is.na(pr)){
        stopifnot(pr>0, pr<1)
        pr.odds <- pr/(1-pr)
        PPV.odds <- pr.odds*LRplus
        PPV <- PPV.odds/(1+PPV.odds)
        NPV.inv.odds <- pr.odds*LRminus
        NPV <- rev(1/(1+NPV.inv.odds))
        res[5, 2:4] <- PPV
        res[6, 2:4] <- NPV
        res[7, 3] <- pr
    }
    if(is.na(pr) & (cohort==TRUE)){
        res[5, 2:4] <- wilson(x=tp, n=tp+fp, conf.level = conf.level)
        res[6, 2:4] <- wilson(x=tn, n=tn+fn, conf.level = conf.level)
        res[7, 2:4] <- wilson(x=tp+fn, n=tp+tn+fp+fn, conf.level = conf.level)
    }

    res[, 1] <- c("Sensitivity", "Specificity", "LRplus", "LRminus", "PPV", "NPV", "prevalence")
    res <- res[,c(1,3,2,4)]
    return(res)
    
}
