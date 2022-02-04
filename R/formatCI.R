################################################################################
### Part of the R package "biostatUZH".
### Free software under the terms of the GNU General Public License (version 2
### or later) a copy of which is available at http://www.R-project.org/Licenses
###
### Copyright (C) 2015 Sebastian Meyer
### Time-stamp: <[formatCI.R] 2015-02-11 18:48 (CET) by SM>
################################################################################


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
