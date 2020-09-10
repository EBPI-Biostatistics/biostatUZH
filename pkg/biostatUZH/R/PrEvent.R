## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function PrEvent
PrEvent <- function(HR, a.length, f.length, dist="exp", lambda=NULL,
                    shape=NULL, survfit.ref=NULL, alloc.ratio=1, method="exact") {
    stopifnot(dist %in% c("exp","exponential","weib","weibull","nonp","non-parametric"))
    stopifnot(HR>0)
    p1 <- alloc.ratio/(alloc.ratio+1)
    p2 <- 1 - p1
    a <- a.length
    f <- f.length
    if(dist %in% c("exp","exponential")){
        stopifnot(lambda>0)
        if(method=="exact"){
            pr.event.cnt <- 1-1/a*integrate(function(x){
                1-pexp(x,lambda)
            }, lower = f, upper = a+f)$value
            pr.event.trt <- 1-1/a*integrate(function(x){
                (1-pexp(x,lambda))^HR
            }, lower = f, upper = a+f)$value
        } else {
            S1 <- 1-pexp(f,lambda)
            S2 <- 1-pexp(.5*a+f,lambda)
            S3 <- 1-pexp(a+f,lambda)
        }
    }
    if(dist %in% c("weib","weibull")){
        stopifnot(!is.null(lambda), !is.null(shape))
        if(method=="exact"){
            pr.event.cnt <- 1-1/a*integrate(function(x){
                1-pweibull(x,shape,lambda)
            }, lower = f, upper = a+f)$value
            pr.event.trt <- 1-1/a*integrate(function(x){
                (1-pweibull(x,shape,lambda))^HR
            }, lower = f, upper = a+f)$value
        }
        else {
            S1 <- 1-pweibull(f,shape,lambda)
            S2 <- 1-pweibull(.5*a+f,shape,lambda)
            S3 <- 1-pweibull(a+f,shape,lambda)
        }
    }
    if(dist %in% c("nonp","non-parametric")){
        stopifnot(!is.null(survfit.ref))
        if(method=="exact") {
            method <- "approx"}
        time <- survfit.ref$time
        surv <- survfit.ref$surv
        index1 <- tail(which(f >= time),1)
        index2 <- tail(which(.5*a+f >= time),1)
        index3 <- tail(which(a+f >= time),1)
        if(is.na(index3)){
            stop("accrual + follow-up length longer than duration of survfit.ref")}
        S1 <- surv[index1]
        S2 <- surv[index2]
        S3 <- surv[index3]
    }
    if(method %in% c("approx","approximate")){
        pr.event.cnt <- 1-1/6*(S1+4*S2+S3)
        pr.event.trt <- 1-1/6*(S1^HR+4*S2^HR+S3^HR)
    }
    pr.event <- p1*pr.event.cnt + p2*pr.event.trt
    return(pr.event)
}
