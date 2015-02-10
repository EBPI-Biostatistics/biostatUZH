formatPercent <- function(x, digits = 1){
  paste(formatC(x * 100, digits = digits, format = "f"), "%", sep = "")
}
