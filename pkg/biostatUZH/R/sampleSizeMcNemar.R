sampleSizeMcNemar <- function(p1, p2, alpha = 0.05, power = 0.8){

## -------------------------------------------------
## sample size for McNemar's test according to
## Lachenbruch (1992), On the sample size for studies based
## upon McNemar's test, Statistics in Medicine, 11, 1521-1525
## -------------------------------------------------

if (p1 > p2){cat("p1 must be smaller than p2")}
if (p2 < 0.5){cat("p2 must be greater than 0.5")}

pp1 <- p1
p1p <- p2
p11 <- sort(seq(min(p1p, pp1), p1p + pp1 - 1, by = -10^-4))
s <- (pp1 - p11) / (pp1 + p1p - 2 * p11)
nl <- ceiling(0.25 * (qnorm(alpha / 2) + qnorm(1 - power)) ^ 2 / (0.5 - s) ^ 2  / (abs(pp1 + p1p - 2 * p11)))
N <- nl[c(1, median(1:length(nl)), length(nl))]
names(N) <- c("N_l min", "N_l mid", "N_l max")

return(N)
}
