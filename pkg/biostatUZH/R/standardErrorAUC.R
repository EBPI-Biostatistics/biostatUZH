
# we follow the second method described in Altman et al, p. 113, 
standardErrorAUC <- function(cases, controls){
    
    ncases <- length(cases)
    ncontrols <- length(controls)
    # non-disease placement values of cases
    C <- rep(NA, ncases)
    # disease placement values of controls
    R <- rep(NA, ncontrols)

    for(i in 1:ncases)
        C[i] <- mean(as.numeric(controls<cases[i])+0.5*as.numeric(controls==cases[i]))
    for(j in 1:ncontrols)
        R[j] <- mean(as.numeric(cases>controls[j])+0.5*as.numeric(cases==controls[j]))
    auc.se <- sqrt((var(R)/ncontrols + var(C)/ncases))

    return(auc.se)
}
