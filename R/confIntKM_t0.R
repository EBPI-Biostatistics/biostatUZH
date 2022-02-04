confIntKM_t0 <- function(time, event, t0, conf.level = 0.95){

alpha <- 1 - conf.level

# compute ci according to the formulas in Huesler and Zimmermann, Chapter 21
obj <- survfit(Surv(time, event) ~ 1, conf.int = 1 - alpha, conf.type = "plain", type = "kaplan", error = "greenwood", conf.lower = "peto")
St <- summary(obj)$surv
t <- summary(obj)$time
n <- summary(obj)$n.risk
res <- matrix(NA, nrow = length(t0), ncol = 3)

for (i in 1:length(t0)){
    ti <- t0[i]
    if (min(t) > ti){res[i, ] <- c(1, NA, NA)}
    if (min(t) <= ti){    
        if (ti %in% t){res[i, ] <- rep(NA, 3)} else {
            Sti <- min(St[t < ti])
            nti <- min(n[t < ti])        
            Var.peto <- Sti ^ 2 * (1 - Sti) / nti
            Cti <- exp(qnorm(1 - alpha / 2) * sqrt(Var.peto) / (Sti ^ (3 / 2) * (1 - Sti)))
            ci.km <- c(Sti / ((1 - Cti) * Sti + Cti), Cti * Sti / ((Cti - 1) * Sti + 1))
            res[i, ] <- c(Sti, ci.km)}
    }
} # end for

res <- cbind(t0, res)
dimnames(res)[[2]] <- c("t0", "S at t0", "lower.ci", "upper.ci")
return(res)
}
