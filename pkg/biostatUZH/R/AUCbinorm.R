AUCbinorm <- function(data, x){
    dat <- data[x, ]
    cases <- dat[dat[, 2] == 1, 1]
    controls <- dat[dat[, 2] == 0, 1]

    a <- (mean(cases) - mean(controls)) / sd(cases)
    b <- sd(controls) / sd(cases)

    auc <- pnorm(a * (b ^ 2 + 1) ^ (-1 / 2))
    return(auc)
}
