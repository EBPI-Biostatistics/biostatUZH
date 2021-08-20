#' Cochran-Mantel-Haenszel Chi-Squared Test for Count Data
#' 
#' Compute Cochran-Mantel-Haenszel chi-squared test of the null that two
#' nominal variables are conditionally independent in each stratum, assuming
#' that there is no three-way interaction.
#' 
#' 
#' @param exposure Binary variable coding whether patient was exposed (1) or
#' not (0).
#' @param outcome Binary variable coding outcome of patient.
#' @param stratum Factor object with at least 2 levels identifying to which
#' stratum the corresponding elements in \code{exposure} and \code{outcome} belong.
#' @return \item{tab}{Table that counts numbers of strata for each case-control
#' combination.}
#' \item{test.stat}{Test statistic for \eqn{chi^2} test.}
#' \item{p.val}{\eqn{p}-value of \eqn{chi^2} test.}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @seealso Similar functionality is provided in \code{\link[stats]{mantelhaen.test}}
#' and \code{\link[survival]{clogit}}. See the examples below for a
#' comparison.
#' @references Agresti, A. (2002). \emph{Categorical data analysis}. Wiley, New York,
#' @keywords htest survival
#' @examples
#' 
#' # generate data
#' set.seed(1977)
#' data <- data.frame(exposure = rep(c(1, 0, 0, 0, 0), 41),
#'                    outcome = sample(x = c(rep(1, 62), rep(0, 5 * 41 - 62))),
#'                    stratum = rep(1:41, each = 5))
#' 
#' # via conditional logistic regression
#' logreg <- clogit(outcome ~ exposure + strata(stratum), method = "approximate",
#'                  data = data)
#' summary(logreg)
#' 
#' # R function in library 'stats'
#' mh <- with(data, mantelhaen.test(x = outcome, y = exposure, z = stratum))
#' 
#' # this function
#' mH <- with(data, mantelHaenszel(exposure = exposure,
#'                                 outcome = outcome,
#'                                 stratum = stratum))
#' 
#' # compare p-values
#' summary(logreg)$coef[5]
#' mh$p.value
#' mH$p.val
#' 
#' @export
mantelHaenszel <- function(exposure, outcome, stratum){
    stopifnot(is.numeric(exposure),
              length(exposure) > 0,
              is.finite(exposure),
              exposure %in% c(0, 1),
              is.numeric(outcome),
              length(outcome) == length(exposure),
              is.finite(outcome),
              outcome %in% c(0, 1))
    stratum <- as.numeric(stratum)
    stopifnot(length(stratum) == length(exposure),
              is.finite(stratum),
              length(unique(stratum)) >= 2)
    
    ss <- unique(stratum)
    tab <- matrix(0, ncol = table(stratum)[1], nrow = 2)
    M <- matrix(ncol = 3, nrow = length(ss))
    
    for (s in ss){
        resp.strat <- outcome[stratum == s]
        expo.strat <- exposure[stratum == s]
        
                                        # generate Mantel-Haenszel table
        resp.case <- resp.strat[expo.strat == 1]
        resp.control <- sum(resp.strat[expo.strat == 0])
        tab[2 - resp.case, resp.control + 1] <- 1 + tab[2 - resp.case, resp.control + 1]
        
                                        # manual computation of MH test statistic
        N <- matrix(NA, 2, 2)
        N[1, 1] <- sum((resp.strat == 1) * (expo.strat == 1))
        N[1, 2] <- sum((resp.strat == 1) * (expo.strat == 0))
        N[2, 1] <- sum((resp.strat == 0) * (expo.strat == 1))    
        N[2, 2] <- sum((resp.strat == 0) * (expo.strat == 0))
        
        m11 <- sum(N[1, ]) * sum(N[, 1]) / sum(N)
        Vn11 <- sum(N[1, ]) * sum(N[2, ]) * sum(N[, 1]) * sum(N[, 2]) / (sum(N) ^ 2 * (sum(N) - 1))           
        
        M[s, ] <- c(N[1, 1], m11, Vn11)
    }
    
    test.stat <- (abs(sum(M[, 1] - M[, 2])) - 0.5) ^ 2 / sum(M[, 3])
    p.val <- 1 - pchisq(test.stat, df = 1)
    
    dimnames(tab) <- list(c("Case exposure yes", "Case exposure no"), rep(0:4))
    list("tab" = tab, "test.stat" = test.stat, "p.val" = p.val)
}







