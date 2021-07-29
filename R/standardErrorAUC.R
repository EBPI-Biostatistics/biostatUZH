# we follow the second method described in Altman et al, p. 113, 


#' Standard Error of AUC
#' 
#' Computes the standard error of the area under the curve.
#' 
#' 
#' @param cases Values of the continuous variable for the cases.
#' @param controls Values of the continuous variable for the controls.
#' @return The standard error.
#' @author Leonhard Held
#' @references The computation follows Chapter 10 in
#' 
#' Altman, D.G., Machin, D., Bryant, T.N. and Gardner, M.J. (2001).
#' \emph{Statistics with confidence}. 2nd Edition, 2000. BMJ Books.
#' @keywords univar htest
#' @export
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
