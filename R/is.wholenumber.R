#' Check if (close to) whole number 
#'
#' Check what elements of a vector are (close to) whole numbers by
#' \code{abs(x - round(x)) < tol}.
#'  
#' @param x Numeric vector
#' @param tol Numeric vector of length 1 specifying the tolerance.
#' Default is \code{.Machine$double.eps^0.5}.
#' @examples
#' is.wholenumber(1L:5L)
#' is.wholenumber(1:5)
#' is.wholenumber(1:5 + .Machine$double.eps^0.5 / 2)
#' is.wholenumber(1:5 + .Machine$double.eps^0.5)
#' @export
is.wholenumber <- function (x, tol = .Machine$double.eps^0.5){
    stopifnot(is.numeric(x), length(x) >= 1,
              is.finite(x),
              is.numeric(tol), length(tol) == 1,
              is.finite(tol), tol > 0)
    abs(x - round(x)) < tol
}
