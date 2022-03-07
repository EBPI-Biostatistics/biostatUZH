#' Format p-Values
#' 
#' Formats p-values according to:
#' \itemize{
#' \item p-values \code{< break.eps} are written as \code{"< break.eps"}.
#' \item p-values \code{>= break.eps} but \code{< break.middle} have 1 digit.
#' \item p-values \code{>= break.middle} have 2 digits.
#' }
#' 
#' @param x Numeric vector of p-values.
#' @param break.eps Numeric vector of length 1.
#' @param break.middle Numeric vector of length 1.
#' @param na.form Character representation of \code{NA}s.
#' @param ... Additional arguments passed to \code{\link{format}}.
#' @return A vector of \code{length(x)} with formatted p-values.
#' @author Sina Rueeger and Sebastian Meyer
#' @seealso the \pkg{base} function \code{\link{format.pval}},
#' \code{\link[Hmisc]{format.pval}} in package \pkg{Hmisc},
#' \code{\link[reporttools]{formatPval}} in package \pkg{reporttools},
#' \code{\link[surveillance]{formatPval}} in package \pkg{surveillance}
#' @examples
#' 
#' x <- c(1e-8, 0.00568, 0.0345, 0.885)
#' biostatUZH::formatPval(x)  # "< 0.0001" "0.006" "0.035" "0.89"
#' 
#' ## compare to formatting of other packages
#' if (requireNamespace("reporttools")) {
#'     reporttools::formatPval(x) # "< 0.0001" "0.0057" "0.03" "0.88"
#' }
#' if (requireNamespace("Hmisc")) {
#'     Hmisc::format.pval(x)  # "0" "0.00568" "0.03450" "0.88500"
#' }
#' if (requireNamespace("surveillance")) {
#'     surveillance::formatPval(x)  # "<0.0001" "0.0057" "0.035" "0.89"
#' }
#' 
#' ## adapt break.middle
#' biostatUZH::formatPval(x, break.middle = 0.001)
#' 
#' @export
formatPval <- function(x, break.eps = 1e-04, break.middle = 0.01, na.form = "NA", ...)
{

    stopifnot(is.numeric(x), length(x) > 0L,
              0 <= x[!is.na(x)], x[!is.na(x)] <= 1,
              is.numeric(break.eps), length(break.eps) == 1,
              is.finite(break.eps), break.eps > 0,
              is.numeric(break.middle), length(break.middle) == 1L,
              is.finite(break.middle), break.middle > 0,
              is.character(na.form), length(na.form) == 1L)

    format1Pval <- function (pv) {
        if (is.na(pv)) {
            na.form 
        } else if (pv < break.eps) {
            paste("<", format(break.eps, scientific = FALSE))
        } else {
            largep <- pv >= break.middle
            format(pv, digits = 1L + largep, nsmall = 1L + largep, scientific = FALSE, ...)
        }
    }
    vapply(X = x, FUN = format1Pval, FUN.VALUE = character(1L), USE.NAMES = TRUE)
}
