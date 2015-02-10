
## =============
## format P-Value
## =============

## x: vector of p-values to format
## Formatting rules:
## * p-values lower than break.eps are formatted as: "< break.eps"
## * p-values lower than break.middle and larger than break.eps have ONE digit
## * p-values larger than break.middle have TWO digits

formatPval <- function(x, break.eps = 1e-04, break.middle = 0.01, ...)
{
    sapply(x, function (xi) {
        if (xi < break.eps) {
            paste("<", format(break.eps, scientific=FALSE))
	} else {
            largep <- xi >= break.middle
            format(xi, digits=1+largep, nsmall=1+largep, scientific=FALSE, ...)
        }
    }, simplify = TRUE, USE.NAMES = TRUE)
}
