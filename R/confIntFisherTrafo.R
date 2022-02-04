confIntFisherTrafo <- function(var1, var2, pp = c(0.025, 0.975), meth = "spearman", type = "t"){

# default values
estimate <- NA
ci1 <- NA
p.val <- NA
p2 <- NA

x <- cbind(var1, var2)
x <- x[is.na(apply(x, 1, sum)) == 0, ]
x <- matrix(x, ncol = 2)
n <- nrow(x)

if (n >= 2){
    corEst <- cor.test(x[, 1], x[, 2], use = "pairwise", method = meth)
    estimate <- corEst$estimate
    p2 <- corEst$p.value
    }

if (n >= 4){
    #cor.est <- cor(x, use = "pairwise", method = meth)[1, 2]
    rho <- 0.5 * log((1 + corEst$estimate) / (1 - corEst$estimate))

    if (type == "z"){
        quant <- qnorm(pp)
        p.val <- pnorm(abs(rho) * sqrt(n - 3))}
    if (type == "t"){
        quant <- qt(pp, df = n - 3)
        p.val <- pt(abs(rho) * sqrt(n - 3), df = n - 3)}

    ci1 <- rho + quant / sqrt(n - 3)
    ci1 <- (exp(2 * ci1) - 1) / (exp(2 * ci1) + 1)

    # does the same
    # ci1 <- tanh(atanh(corEst$estimate) + quant / sqrt(n - 3))

    p.val <- 2 * (1 - p.val)
    }

res <- list("estimate" = as.numeric(estimate), "ci" = ci1, "p.value" = as.numeric(p.val), "n" = n, "p2" = p2)
return(res)
}








