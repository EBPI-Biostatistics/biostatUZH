#' Test for qualitative interaction by Gail and Simon
#' 
#' @description This function computes the p-value of the Gail-Simon test for
#' qualitative interaction from effect estimates and the corresponding 
#' standard errors in a number of subgroups
#'
#' @param thetahat A numeric vector of length > 1. Denotes the effect estimates
#' in the different subgroups.
#' @param se A numeric vector of the same length as \code{thetahat}. Contains 
#' the standard errors of the corresponding effect estimates in each subgroup.
#'
#' @return A numeric vector of length 1 containing the p-value according 
#' to the Gail-Simon test.
#' 
#' @references 
#' Gail, M., & Simon, R. (1985). Testing for Qualitative Interactions between 
#' Treatment Effects and Patient Subsets. \emph{Biometrics, 41}(2), 361–372. 
#' \url{https://doi.org/10.2307/2530862}
#' 
#' @author Leonhard Held
#' 
#' @export
#'
#' @examples
#' rd <- c(0.163, -0.114, -0.047, -0.151)
#' rd.se <- c(0.0788, 0.0689, 0.0614, 0.0547)
#' ## Gail and Simon test
#' gailSimon(thetahat = rd, se = rd.se)
gailSimon <- function(thetahat, se){
  stopifnot(is.numeric(thetahat), length(thetahat) > 1,
            !any(is.na(thetahat)), 
            is.numeric(se), length(se) > 1,
            !any(is.na(se)),
            length(thetahat) == length(se))
  nSubgroups <- length(thetahat)
  myrange <- 1L:(nSubgroups - 1L)
  Qplus <- sum((thetahat >= 0) * (thetahat/se)^2)
  Qminus <- sum((thetahat < 0) * (thetahat/se)^2)
  minQ <- min(Qplus, Qminus)
  pval <- sum(stats::dbinom(x = myrange, size = nSubgroups - 1L, prob = 0.5) * (1 - stats::pchisq(minQ, df = myrange)))
  return(pval)
}
