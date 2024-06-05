#' Square-and-Add Method for Confidence Interval of Difference in Point Estimates
#'
#' This function computes the confidence interval for the difference between 
#' two point estimates using the square-and-add method.
#'
#' @param theta1 Numeric value of the point estimate for group 1.
#' @param lower1 Numeric value of the lower bound of the confidence interval for group 1.
#' @param upper1 Numeric value of the upper bound of the confidence interval for group 1.
#' @param theta2 Numeric value of the point estimate for group 2.
#' @param lower2 Numeric value of the lower bound of the confidence interval for group 2.
#' @param upper2 Numeric value of the upper bound of the confidence interval for group 2.
#' 
#' @return A list with the entries: \itemize{ \item difference: Estimated difference (group 1 - group 2).
#' \item CI: Lower and upper bound of the confidence interval for the difference.}
#' 
#' @details
#' The function returns a confidence interval of the same level 
#' as the confidence intervals for theta1 and theta2
#' 
#' @author Manuel Pfister, Leonhard Held, Charlotte Micheloud
#' 
#' @references 
#' Newcombe, R.G. (1998). Interval estimation for the difference 
#' between independent proportions: comparison of eleven methods, 
#' \emph{Statistics in Medicine}, \bold{17}, 873--890.
#' 
#' Newcombe, R.G. (2013). \emph{Confidence Intervals for Proportions 
#' and Related Measures of Effect Size}, CRC Press, 
#' 2nd Edition, Chapter 7.3.
#' 
#' @examples
#' ci1 <- wilson(x=11, n=12)
#' ci2 <- wilson(x=1, n=12)
#' confIntSquareAdd(theta1 = ci1[2], lower1=ci1[1], upper1=ci1[3],
#'                  theta2 = ci2[2], lower2=ci2[1], upper2=ci2[3]) 
#' @export

confIntSquareAdd <- function(theta1, lower1, upper1, theta2, lower2, upper2) {
  # Error handling
  if (!is.numeric(theta1) || !is.numeric(lower1) || !is.numeric(upper1) ||
      !is.numeric(theta2) || !is.numeric(lower2) || !is.numeric(upper2)) {
    stop("All inputs must be numeric values.")
  }

  if (lower1 > theta1 || upper1 < theta1 || lower2 > theta2 || upper2 < theta2) {
    stop("Lower bound must be less than or equal to the point estimate, 
         and upper bound must be greater than or equal to the point estimate.")
  }

  # Compute the difference in point estimates
  diff <- theta1 - theta2

  # Compute the confidence interval using the square-and-add method
  diff.lower <- diff - sqrt((theta1 - lower1)^2 + (upper2 - theta2)^2)
  diff.upper <- diff + sqrt((upper1 - theta1)^2 + (theta2 - lower2)^2)

  res <- data.frame(matrix(c(diff.lower, diff.upper), ncol = 2))
  names(res) <- c("lower", "upper")
  
  list("difference" = diff, "CI" = res)
}
