## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function NumEvents
NumEvents <- function(HR, sig.level=0.05, power=NULL, n.events=NULL,
                     alloc.ratio=1, non.inf.margin=NULL, type="sup",
                     alternative="two.sided"){
    stopifnot(HR>0)
    logHR <- log(HR)
    stopifnot(alternative %in% c("two.sided","one.sided"))
    tails <- ifelse(alternative=="one.sided",1,2)
    z_alpha <- qnorm(sig.level/tails)
    p <- alloc.ratio/(alloc.ratio+1)
    p.factor <- p*(1-p)
    ## calculate # events
    if(is.null(n.events)){
        z_beta <- qnorm(1-power)
        c <- (z_alpha+z_beta)^2
        ## superiority case
        if(type %in% c("superiority","sup")){
            d <- c/(p.factor*logHR^2)
            return(ceiling(d))
        }
        ## non-inferiority case
        if(type %in% c("non-inferiority","noninf")){
            stopifnot(!is.null(non.inf.margin))
            d <- c/(p.factor*(logHR - non.inf.margin)^2)
            return(ceiling(d))
        }
    }
    ## calculate power
    if(!is.null(n.events)){
        ## superiority case
        if(type %in% c("superiority","sup")){
            power <- 1 - pnorm(logHR*sqrt(n.events*p.factor)-z_alpha)
            return(power)
        }
        ## non-inferiority case
        if(type %in% c("non-inferiority","noninf")){
            stopifnot(!is.null(non.inf.margin))
            power <- 1 - pnorm((logHR-non.inf.margin)*sqrt(n.events*p.factor)-z_alpha)
            return(power)
        }
    }
    stop("Wrong distribution")
}
