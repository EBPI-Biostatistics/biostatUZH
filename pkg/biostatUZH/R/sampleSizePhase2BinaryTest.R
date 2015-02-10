sampleSizePhase2BinaryTest <- function(FPF0 = 0.2, TPF0 = 0.75, FPF1 = 0.05, TPF1 = 0.9, alpha = 0.05, power = 0.9){

    beta <- 1 - power
    astar <- 1 - sqrt(1 - alpha)
    bstar <- 1 - sqrt(1 - beta)

    za <- qnorm(1 - astar)
    zb <- qnorm(1 - bstar)

    nd <- (za * sqrt(TPF0 * (1 - TPF0)) + zb * sqrt(TPF1 * (1 - TPF1))) ^ 2 / (TPF1 - TPF0) ^ 2
    nd2 <- (za * sqrt(FPF0 * (1 - FPF0)) + zb * sqrt(FPF1 * (1 - FPF1))) ^ 2 / (FPF1 - FPF0) ^ 2

    res <- list("n.cases" = ceiling(nd), n.controls = ceiling(nd2))
    return(res)
}
