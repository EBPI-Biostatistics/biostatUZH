#' Power calculations for two sample z tests
#'
#' @param n Number of observations (per group)
#' @param delta True difference in means
#' @param sd Standard deviation
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param type String specifying the type of t test. Can be abbreviated.
#' @param alternative One- or two-sided test. Can be abbreviated.
#'
#' @details Exactly one of the parameters \code{n}, \code{delta}, 
#' \code{power}, \code{sd}, and \code{sig.level} must be passed as \code{NULL}, 
#' and that parameter is determined from the others. Notice that the last two 
#' have non-\code{NULL} defaults, so \code{NULL} must be explicitly passed if 
#' you want to compute them.
#' 
#' @return Object of class "power.htest", a list of the arguments 
#' (including the computed one) augmented with method and note elements.
#' 
#' @author Felix Hofmann
#' @export
#'
#' @examples
#' 
#' # Calculate sample size
#' delta = 0.25
#' sd = 0.4
#' sig.level = 0.01
#' power = 0.95
#' power.z.test(delta = delta, sd = sd, sig.level = sig.level, power = power)
#' 
#' # Calculate the effect size
#' n = 92
#' power.z.test(power = power, sd = sd, sig.level = sig.level, n = n)
#' 
#' # Calculate the standard deviation
#' power.z.test(power = power, delta = delta, sig.level = sig.level, n = n,
#'              sd = NULL)
#' 
#' # Calculate the type I error
#' power.z.test(power = power, delta = delta, sig.level = NULL, n = n,
#'              sd = sd)
#' 
#' 
#' # Calculate power
#' power.z.test(delta = delta, sd = sd, sig.level = sig.level, n = n)
#' 
power.z.test <- function (n = NULL, 
                          delta = NULL, 
                          sd = 1, 
                          sig.level = 0.05, 
                          power = NULL,
                          type = c("two.sample", "one.sample"),
                          alternative = c("two.sided", "one.sided")
                          ){
  
  # Check exactly one of n, delta, sd, power, sig.level is NULL
  # Adapted from stats::power.t.test
  if (sum(vapply(list(n, delta, sd, power, sig.level), is.null, 
                 logical(1L))) != 1)
    stop("exactly one of 'n', 'delta', 'sd', 'power', and 'sig.level' must be NULL")
  # Check that n is larger than 1
  if(!is.null(n) && n < 1)
    stop("n must be larger than 1")
  # Check that sd is positive
  if(!is.null(sd) && sd <= 0)
    stop("sd must be positive")
  
  # Check that sig.level and power are either NULL or between 0 and 1
  assert_NULL_or_prob(sig.level)
  assert_NULL_or_prob(power)
  
  # set alternative and type for now
  type <- match.arg(type)
  tsample <- switch(type, one.sample = 1, two.sample = 2)# , paired = 1)
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  
  if(is.null(n)){
    n <- tsample * ((qnorm(power) + qnorm(1 - sig.level/tside)) * sd / delta)^2
  } else if(is.null(delta)){
    delta <- sd * sqrt(tsample / n) * (qnorm(power) + qnorm(1 - sig.level/tside))
  } else if(is.null(sd)){
    sd <- abs(delta) * sqrt(n) / (sqrt(tsample) * (qnorm(power) + qnorm(1 - sig.level/tside)))
  } else if(is.null(sig.level)){
    sig.level <- tside * (1 - pnorm(sqrt(n/tsample) * abs(delta)/sd - qnorm(power)))
  } else if(is.null(power)){
    power <- pnorm(sqrt(n/tsample) * abs(delta) / sd - qnorm(1 - sig.level/tside))
  }
  
  # Create notes (copied from stats::power.t.test)
  NOTE <- switch(type, 
                 paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs", 
                 two.sample = "n is number in *each* group", 
                 NULL)
  METHOD <- paste(switch(type, 
                         one.sample = "One-sample", 
                         two.sample = "Two-sample", 
                         paired = "Paired"), 
                  "z test power calculation")
  
  # Return object (copied from stats::power.t.test)
  structure(list(n = n, delta = delta, sd = sd, sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")

}

# The following function is a copy of stats:::assert_NULL_or_prob
# We use a copy of this function because we do not want to import internal 
# functions from other packages

#' @noRd 
assert_NULL_or_prob <- function(x){
  if(!is.null(x) && (!is.numeric(x) || any(0 > x | x > 1)))
    stop(gettextf("'%s' must be numeric in [0, 1]", deparse1(substitute(x))), 
         call. = FALSE)
  invisible(x)
}

