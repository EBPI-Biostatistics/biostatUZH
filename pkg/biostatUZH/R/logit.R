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
