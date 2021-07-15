logit <- function (x)
{
    if (any(omit <- is.na(x) | x <= 0 | x >= 1)) {
        is.na(x) <- omit
        if (any(!omit))
            x[!omit] <- Recall(x[!omit])
        x
    } else qlogis(x) # = log(x/(1 - x))
}

ilogit <- function (x) plogis(x) # = exp(x)/(1 + exp(x))
expit <- ilogit



logit_fg <- function(x)
{
    stopifnot(!is.null(x))
    ind <- x <= 0 | x >= 1
    if(any(ind)){
        warning("Values not in [0,1] set to NA")
        x[ind] <- NA
    }
    qlogis(x) # = log(x/(1 - x))
}



x <- c(1,2,.5,NA, NULL)
all.equal(logit(x), logit_fg(x))
