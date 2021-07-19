#' Logit and inverse-logit function
#' 
#' Compute logit and inverse logit function.
#' 
#' Used by \code{\link{faganPlot}}, \code{\link{faganLine}}.
#' 
#' 
#' @aliases logit ilogit expit
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
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
