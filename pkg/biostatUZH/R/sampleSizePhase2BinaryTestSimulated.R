sampleSizePhase2BinaryTestSimulated <- function(FPF0, TPF0, FPF1, TPF1, alpha = 0.05, power = 0.9, M = 5000, range = 15, print = TRUE, type = c("wilson", "exact")[2]){

## simulate power for a Phase 2 diagnostic test using rectangular confidence region
## based on one-sided exact confidence limits
## Example 8.1 in Pepe (2001), The Statistical Evaluation Medical Tests for Classification and Prediction
## initial guess from normal approximation

require(biostatUZH)   # to provide sampleSizePhase2BinaryTest()

## significance level for rectangular region
beta <- 1 - power
alphastar <- 1 - (1 - alpha) ^ (1 / 2)

## sample sizes from normal approximation
wald <- sampleSizePhase2BinaryTest(FPF0, TPF0, FPF1, TPF1, alpha, power)

## vector of sample sizes to simulate power for
nDs <- wald$n.cases + (-range):range
nDbars <- wald$n.controls + (-range):range

## result matrix
powers <- matrix(NA, ncol = 3, nrow = length(nDs))
colnames(powers) <- c("nD", "nDbars", "power")

## define confidence interval to be used
if (identical(type, "wilson")){confint <- wilson}
if (identical(type, "exact")){confint <- clopperPearson}

for (s in 1:length(nDs)){

    nD <- nDs[s]
    nDbar <- nDbars[s]

    ## simulate number of diseased cases and controls
    nDplus <- rbinom(M, nD, TPF1)
    nDplusbar <- rbinom(M, nDbar, FPF1)

    decisions <- rep(NA, M)

    for (i in 1:M){

        ## estimate TPF and FPF for each simulated number
        TPFhat <- nDplus[i] / nD
        FPFhat <- nDplusbar[i] / nDbar

        ## compute exact interval for each proportion
        ciTPF <- confint(nDplus[i], nD, conf.level = 1 - 2 * alphastar)
        ciFPF <- confint(nDplusbar[i], nDbar, conf.level = 1 - 2 * alphastar)

        ## decide on H0 or H1 using a rectangular confidence region using *one-sided* exact confidence limits
        ind1 <- ciTPF[1] >= TPF0
        ind2 <- ciFPF[3] <= FPF0
        decisions[i] <- 0
        if (ind1 & ind2){decisions[i] <- 1}
    }

    powers[s, ] <- c(nD, nDbar, sum(decisions) / M)
    if (identical(print, TRUE)){print(paste("s = ", s, " of ", length(nDs), " done", sep = ""))}
}

firstReallyAbovePower <- max((1:nrow(powers))[powers[, "power"] < 1 - beta]) + 1
simul <- powers[firstReallyAbovePower, ]

res <- list("wald.n.cases" = wald$n.cases, "wald.n.controls" = wald$n.controls, "simul.n.cases" = as.numeric(simul["nD"]), "simul.n.controls" = as.numeric(simul["nDbars"]),
    "simul.power" = as.numeric(simul["power"]), "alphastar" = alphastar, "simul.powers" = powers)
return(res)
}
