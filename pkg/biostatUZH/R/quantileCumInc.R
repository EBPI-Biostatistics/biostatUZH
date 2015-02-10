quantileCumInc <- function(time, event, group, quant = 0.5)
{
if (!requireNamespace("cmprsk")) stop("requires cmprsk::cuminc()")

t.low0 <- unlist(lapply(split(time, interaction(group, event, sep = " ")), min))
t.up0 <- unlist(lapply(split(time, interaction(group, event, sep = " ")), max))

est0 <- cmprsk::timepoints(cmprsk::cuminc(time, event, group = group, rho = 0, cencode = 0), times = mean(time))$est
res <- matrix(NA, nrow = nrow(est0), ncol = 1)
rownames(res) <- rownames(est0)
colnames(res) <- paste(quant, "-quantile", sep = "")

for (r in 1:nrow(est0)){
    tmp <- function(x, time, event, group){cmprsk::timepoints(cmprsk::cuminc(time, event, group = group, rho = 0, cencode = 0), times = x)$est[r, 1] - quant}
        
    low <- tmp(x = t.low0[r], time, event, group)
    up <- tmp(x = t.up0[r], time, event, group)

    i2 <- 0
    i4 <- 0
    i5 <- 0

    i1 <- is.na(low) == FALSE
    if (i1 == 1){i2 <- low <= 0}
    i3 <- is.na(up) == FALSE
    if (i3 == 1){i4 <- 0 <= up}
    if (i1 * i3 == 1){i5 <- low < up}
    if (i1 * i2 * i3 * i4 * i5 == 1){res[r, 1] <- uniroot(tmp, interval = c(t.low0[r], t.up0[r]), time, event, group)$root}
}

return(res)
}
