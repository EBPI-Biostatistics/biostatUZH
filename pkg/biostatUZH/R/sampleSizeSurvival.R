## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function sampleSizeSurvival
sampleSizeSurvival <- function(HR, a.length, f.length, sig.level=0.05, power=NULL,
                               n = NULL, n.events=NULL, alloc.ratio=1, drop.rate=0,
                               non.inf.margin=NULL, type="sup", dist="exp",
                               lambda=NULL, shape=NULL, survfit.ref=NULL,
                               alternative="two.sided", method="exact") {
    ## probability of event
    pr.event <- PrEvent(HR = HR, a.length = a.length, f.length = f.length,
                        dist = dist, lambda = lambda, shape = shape, method = method,
                        survfit.ref = survfit.ref, alloc.ratio = alloc.ratio)
    ## number of events or power needed (depending on given arguments)
    if(is.null(power) & is.null(n.events) & is.null(n)){
        stop(paste("either the power or the number of events",
                   "or the sample size need to be specified", sep = "\n"))
    }
    ## power calculation
    if(is.null(power)){
        if(is.null(n.events)){
            n.events <- ceiling(pr.event * n)
        }
        power <- NumEvents(HR = HR, sig.level = sig.level, power = power,
                          n.events = n.events, alloc.ratio = alloc.ratio,
                          non.inf.margin = non.inf.margin, type = type,
                          alternative = alternative)
    }
    ## sample size calculation
    if(is.null(n)){
        if(is.null(n.events)){
            n.events <- NumEvents(HR = HR, sig.level = sig.level, power = power,
                                 n.events = n.events, alloc.ratio = alloc.ratio,
                                 non.inf.margin = non.inf.margin, type = type,
                                 alternative = alternative)
        }
        n <- n.events/pr.event
        n <- ceiling(n/(1-drop.rate))}
    ## number of events calculation
    ## (only happens when power and sample given, n.event=NULL)
    if(is.null(n.events)){
        n.events <- NumEvents(HR = HR, sig.level = sig.level, power = power,
                             n.events = n.events, alloc.ratio = alloc.ratio,
                             non.inf.margin = non.inf.margin, type = type,
                             alternative = alternative)
    }
    ## return statement
    str <- structure(list(n=n, HR = HR, power=power, sig.level=sig.level,
                          alternative=alternative, distribution = dist,
                          PrEvent = pr.event, Events=n.events,
                          NOTE="n is the total sample size"))
    return(str)
}
