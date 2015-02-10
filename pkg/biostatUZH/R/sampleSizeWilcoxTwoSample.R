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

