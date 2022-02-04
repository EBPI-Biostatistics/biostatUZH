
## computes the standard error of the difference in AUC of two paired tests
## cases is a m x 2 matrix
## controls is a n x 2 matrix
## each row corresponds to measurements from one case/control

standardErrorAUCDiff <- function(cases, controls){

    ncases <- nrow(cases)
    ncontrols <- nrow(controls)
    # non-disease placement values of cases
    C <- matrix(NA,nrow=ncases, ncol=2)
    # disease placement values of controls
    R <- matrix(NA,nrow=ncontrols, ncol=2)

    for(k in 1:2){
        for(i in 1:ncases)
            C[i,k] <- mean(as.numeric(controls[,k]<cases[i,k])+0.5*as.numeric(controls[,k]==cases[i,k]))
        for(j in 1:ncontrols)
            R[j,k] <- mean(as.numeric(cases[,k]>controls[j,k])+0.5*as.numeric(cases[,k]==controls[j,k]))
    }
    auc.diff.se <- sqrt((var(R[,1]-R[,2])/ncontrols + var(C[,1]-C[,2])/ncases))
    return(auc.diff.se)
}
