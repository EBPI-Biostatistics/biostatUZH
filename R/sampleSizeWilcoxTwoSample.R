#' Compute sample size for two sample Wilcoxon (Mann-Whitney) test
#' 
#' Compute sample size to test the hypothesis that two samples come from the
#' same population against that \eqn{Y}'s tend to be larger than \eqn{X}'s.
#' 
#' Given two independent samples \eqn{X_1, ..., X_m} and \eqn{Y_1, ..., Y_n},
#' we want to test the hypothesis that the two samples come from the same
#' population against that \eqn{Y}'s tend to be larger than \eqn{X}'s.
#' 
#' @param a Significance level of test.
#' @param b Desired power of test.
#' @param c Proportion of observations in group 1: \eqn{c = m / (m + n)}
#' (\eqn{c = 0.5} means equally sized groups).
#' @param pxy A value for the probability \eqn{P(Y > X)}.
#' @param two.sided If \code{TRUE} a two-sided test is assumed, otherwise
#' one-sided.
#' @return \item{m}{Sample size of the first group.} \item{n}{Sample size of
#' the second group.}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Nother, G.E. (1987). Sample Size Determination for Some Common
#' Nonparametric Tests \emph{JASA}, \bold{82}, 644--647.
#' @keywords htest
#' @examples
#' 
#' # compute sample size for some pxy's
#' pxys <- c(0.65, 0.7, 0.75)
#' dat1 <- matrix(ncol = 3, nrow = length(pxys))
#' colnames(dat1) <- c("P(Y > X)", "m", "n")
#' for (j in 1:length(pxys)){dat1[j, ] <- c(pxys[j], 
#'     sampleSizeWilcoxTwoSample(a = 0.05, b = 0.1, c = 0.5, 
#'     pxy = pxys[j], two.sided = TRUE))}
#' dat1
#' 
sampleSizeWilcoxTwoSample <- function(a = 0.05, b = 0.2, c = 0.5, pxy = 0.75, two.sided = TRUE){

# Given two independent samples X_1, ..., X_m and Y_1, ..., Y_n, 
# we want to test the hypothesis that the two samples come from
# the same population against that Y's tend to be larger than X's
# Taken from Nother (1987), JASA 82, 644 - 647
#
# Input:
# a:            significance level
# b:            desired power
# c:            proportion of observations in group 1: c = m / (m + n)
#               (c = 0.5 means equally sized groups)
# pxy:          P(Y > X)  
# two.sided:    if = T --> two.sided, else one sided test 
#
# Output:       vector E N^2 containing
# m:            number of observations in Sample 1
# n:            number of observations in Sample 2
#
# Kaspar Rufibach, December 2007
#
    side <- 2
    if (two.sided == FALSE){side <- 1}
    za <- qnorm(1 - a / side) 
    zb <- qnorm(1 - b)
    N <- (za + zb) ^ 2 / (12 * c * (1 - c) * (pxy - 1 / 2) ^ 2)    
    m <- ceiling(c * N)
    n <- ceiling(N * (1 - c))
    
    return(c("m = " = m, "n = " = n))
}

