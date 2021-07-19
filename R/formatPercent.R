#' Format a Numeric Proportion
#' 
#' Takes a number and formats it as a percentage.  Used by
#' \code{\link{confIntIndependentProportion}}.
#' 
#' 
#' @param x numeric vector of proportions
#' @param digits number of digits
#' @author Leonhard Held
#' @examples
#' 
#' formatPercent(c(0.115, 0.5))  # "11.5%" "50.0%"
#' 
formatPercent <- function(x, digits = 1){
  paste(formatC(x * 100, digits = digits, format = "f"), "%", sep = "")
}
