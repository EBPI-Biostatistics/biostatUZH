## computes the standard error of the difference in AUC of two paired tests
## cases is a m x 2 matrix
## controls is a n x 2 matrix
## each row corresponds to measurements from one case/control



#' Standard Error of a AUC difference
#' 
#' Computes the standard error of the difference in area under the curve for
#' paired samples.
#' 
#' 
#' @param cases Matrix with values of the continuous variable for the cases.
#' First column gives values of test A, second gives values of test B.
#' @param controls Matrix with values of the continuous variable for the
#' controls. First column gives values of test A, second gives values of test
#' B.
#' @return The standard error.
#' @author Leonhard Held
#' @references Pepe, M.S. (2003) \emph{The statistical evaluation of medical
#' tests for classification and prediction}. Oxford: Oxford University Press.
#' @keywords univar htest
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
