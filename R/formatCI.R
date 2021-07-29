################################################################################
### Part of the R package "biostatUZH".
### Free software under the terms of the GNU General Public License (version 2
### or later) a copy of which is available at http://www.R-project.org/Licenses
###
### Copyright (C) 2015 Sebastian Meyer
### Time-stamp: <[formatCI.R] 2015-02-11 18:48 (CET) by SM>
################################################################################




#' Format (Confidence) Intervals
#' 
#' The function \code{formatCI} formats (confidence) intervals with limits
#' either represented as text or math symbols and optionally including a unit.
#' It is a re-implementation of \code{\link[reporttools]{displayCI}} from the
#' \pkg{reporttools} package, allowing for customization of the \code{text}
#' format as well as element-wise \code{digits} and \code{unit} arguments.
#' 
#' 
#' @param x either a numeric vector of length 2 representing a single interval,
#' or a 2-column matrix of several intervals.
#' @param digits the number of digits after the decimal point.
#' @param unit a character string denoting a measurement unit.
#' @param text either a character string referring to one of the predefined
#' text styles (\code{"none"}, \code{"german"}, \code{"english"}), or a
#' character vector of length 3 defining the strings to put before, between,
#' and after the limits.
#' @return A character vector of formatted (confidence) intervals.
#' @author Sebastian Meyer
#' @keywords print character
#' @examples
#' 
#' ## single interval, default style
#' formatCI(c(0.9876, 1.2345))
#' ## "[0.99, 1.23]"
#' 
#' ## matrix of several intervals
#' cis <- matrix(c(
#'     0.0123, 0.456,
#'     -10.003, 5.3,
#'     pi, 2*pi
#' ), nrow = 3, ncol = 2, byrow = TRUE)
#' 
#' ## use only 1 digit
#' formatCI(cis, digits = 1)
#' ## "[0.0, 0.5]" "[-10.0, 5.3]" "[3.1, 6.3]"
#' 
#' ## with customized numbers of digits and units
#' formatCI(cis, digits = 1:3, unit = c("cm", "mm", "m"))
#' ## "[0.0cm, 0.5cm]" "[-10.00mm, 5.30mm]" "[3.142m, 6.283m]"
#' 
#' ## with a textual representation of the limits
#' formatCI(cis, text = "english")
#' ## "from 0.01 to 0.46" "from -10.00 to 5.30" "from 3.14 to 6.28"
#' formatCI(cis, text = c("", " -- ", ""))
#' ## "0.01 -- 0.46" "-10.00 -- 5.30" "3.14 -- 6.28"
#' 
#' @export
formatCI <- function (x, digits = 2, unit = "", text = "none")
{
    ## parse arguments
    x <- if (is.vector(x) && length(x) == 2L) t(x) else as.matrix(x)
    stopifnot(is.numeric(x), ncol(x) == 2, is.vector(text, mode = "character"))
    if (length(text) == 1L) {
        text <- switch(
            text,
            "none" = c("[", ", ", "]"),
            "german" = c("von ", " bis ", ""),
            "english" = c("from ", " to ", ""),
            stop("there is no predefined 'text' style \"", text, "\"")
        )
    } else if (length(text) != 3) {
        stop("'text' must be a single character string or of length 3")
    }
    
    ## format the confidence interval(s)
    fmtlimit <- paste0("%.", digits, "f")
    fmt <- paste0(text[1L], fmtlimit, "%s", text[2L], fmtlimit, "%s", text[3L])
    sprintf(fmt, x[,1L], unit, x[,2L], unit)
}
