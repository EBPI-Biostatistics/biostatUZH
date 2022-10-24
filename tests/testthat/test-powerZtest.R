# Function that tests whether the results are consistent given a concrete
# alternative and type
test_scenarios <- function(delta, sd, sig.level, power, type, alternative,
                           missing){
  
  # calculate n
  args <- list(power.z.test, n = NULL, delta = delta, sd = sd, 
               sig.level = sig.level, power = power, type = type, 
               alternative = alternative)
  res <- eval(as.call(args))
  
  # add n to the argument list
  args$n <- res$n
  
  # test whether the function returns the same object, independent of which
  # parameter is calculated
  comp <- vapply(missing, function(x, args, res){
    args[x] <- list(NULL)                    # set missing parameter to NULL
    out <- eval(as.call(args))               # Calculate it
    isTRUE(all.equal(out, res))              # Check whether it is the same
  }, args = args, res = res, logical(1L))
  
  comp
}

# For each of the scenarios, calculate n, then estimate all the others
# This ensures that estimation gives the same results independent of which
# parameter is estimated.

test_that("power.z.test results are consistent", {
  # Set test values
  delta <- 0.25
  sd <- 0.4
  power <- 0.95
  sig.level <- 0.01
  # Set missing parameters (all except for n)
  missing <- c("delta", "sd", "sig.level", "power")
  # Set the possible scenario combinations
  grid <- expand.grid(
    type = c("one.sample", "two.sample"),               # one or two sample test
    alternative = c("one.sided", "two.sided"),          # one or two sided
    stringsAsFactors = FALSE
  )
  # Loop over all combinations of type and alternative and test consistency
  out <- vapply(seq_len(nrow(grid)), function(x){
    # Set type and alternative
    type <- grid$type[x]
    alternative <- grid$alternative[x]
    # Check whether estimates are consistent
    all(test_scenarios(delta = delta, sd = sd, sig.level = sig.level,
                       power = power, type = type, alternative = alternative,
                       missing = missing))
  }, logical(1L))
  expect_true(all(out))
})
