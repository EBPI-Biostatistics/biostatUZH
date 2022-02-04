mantelHaenszel <- function(exposure, outcome, stratum){

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
res <- list("tab" = tab, "test.stat" = test.stat, "p.val" = p.val)
return(res)
}







