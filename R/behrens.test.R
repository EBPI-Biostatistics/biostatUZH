#' Behrens version of the t-test for unequal variances
#'
#' Behrens version of the t-test for unequal variances.
#' @details Computes Behrens version of the t-test for unequal variances
#' based on approximate solution as in Box & Tiao, 1973, Section 2.5.3;
#' see also Armitage, Berry, Matthews, 2002, Section 4.3.
#' @export
behrens.test <- function(x, ...){
  UseMethod("behrens.test")
}


#' @param x Numeric vector of data values.
#' @param y Numeric vector of data values.
#' @param alternative Character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater", and "less".
#' @param mu A number indicating the difference in means under the
#' null-hypothesis.
#' @param conf.level Confidence level of the interval, default is 0.95.
#' @examples
#' weightWomen <- c(38.9, 61.2, 73.3, 21.8, 63.4, 64.6, 48.4, 48.8, 48.5)
#' weightMen <-   c(67.8, 60  , 63.4, 76,   89.4, 73.3, 67.3, 61.3)
#' 
#' # using the default method
#' behrens.test(x = weightWomen, y = weightMen)
#' behrens.test(x = weightWomen, y = weightMen, alternative = "less")
#' 
#' @rdname behrens.test
#' @exportS3Method behrens.test default
behrens.test.default <- function(x, y,
                                 alternative = c("two.sided", "greater", "less"),
                                 mu = 0,
                                 conf.level = 0.95, ...){
  stopifnot(is.numeric(x), length(x) >= 1L, is.finite(x),
            is.numeric(y), length(y) >= 1L, is.finite(y),
            is.numeric(conf.level), length(conf.level) == 1L,
            is.numeric(mu), length(mu) == 1L,
            0 < conf.level, conf.level < 1,
            !is.null(alternative))
  alternative <- match.arg(alternative)
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  alpha <- 1 - conf.level
  m1 <- mean(x)
  s1 <- stats::sd(x)
  n1 <- sum(!is.na(x))
  m2 <- mean(y)
  s2 <- stats::sd(y)
  n2 <- sum(!is.na(y))
  
  nu1 <- n1 - 1
  nu2 <- n2 - 1
  diff <- m1 - m2
  diff.se <- sqrt(s1^2 / n1 + s2^2 / n2)
  
  cos2.phi <- (s2^2 / n2) / diff.se^2
  sin2.phi <- (s1^2 / n1) / diff.se^2
  
  f1 <- nu2 / (nu2 - 2) * cos2.phi + nu1 / (nu1 - 2) * sin2.phi 
  f2 <- nu2^2 / ((nu2 - 2)^2 * (nu2 - 4)) * cos2.phi^2 + nu1^2 /
    ((nu1 - 2)^2 * (nu1 - 4)) * sin2.phi^2 
  df <- 4 + f1^2 / f2
  a <- sqrt((df - 2) / df * f1)
  diff.se <- a * diff.se
  t.value <- (diff - mu) / diff.se
  
  if(alternative == "two.sided"){
    c <- stats::qt(alpha / 2, df = df, lower.tail = FALSE)
    lower <- diff - diff.se * c
    upper <- diff + diff.se * c
    p <- 2 * stats::pt(abs(t.value), df = df, lower.tail = FALSE)
  } else if(alternative == "less"){
    c <- stats::qt(alpha, df = df, lower.tail = FALSE)
    lower <- -Inf
    upper <- diff + diff.se * c
    p <-  stats::pt(t.value, df = df, lower.tail = TRUE)
  } else {
    c <- stats::qt(alpha, df = df, lower.tail = FALSE)
    lower <- diff - diff.se * c
    upper <- Inf
    p <-  stats::pt(t.value, df = df, lower.tail = FALSE)
  }
  names(t.value) <- "t"
  names(df) <- "df"
  names(mu) <- "difference in means"
  estimate <- c(m1, m2)
  names(estimate) <- c("mean of x", "mean of y")
  cint <- c(lower, upper)
  attr(cint, "conf.level") <- conf.level
  method <- "Behrens' t-test"
  rval <- list(statistic = t.value, parameter = df, p.value = p, 
               conf.int = cint, estimate = estimate, null.value = mu, 
               stderr = diff.se, alternative = alternative, method = method, 
               data.name = dname)
  class(rval) <- "htest"
  rval
}

#' @param formula A formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric 
#' variable giving the data values and \code{rhs} is a factor with two levels giving 
#' the corresponding groups.
#' @param data An optional matrix or data frame 
#' (or similar: see \link[stats]{model.frame}) containing the variables in the 
#' formula formula. By default the variables are taken from \code{environment(formula)}.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param na.action A function which indicates what should happen when the data 
#' contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param ... Further arguments to be passed to or from methods.
#' @rdname behrens.test
#' @examples
#' # using the formula interface
#' weight <- matrix(c(weightWomen, weightMen, rep(1, 9), rep(2, 8)), ncol = 2)
#' colnames(weight) <- c("measure", "sex") 
#' behrens.test(measure ~ sex, data = weight)
#' @exportS3Method behrens.test formula
behrens.test.formula <- function (formula, data, subset, na.action, ...){
  if (missing(formula) || (length(formula) != 3L)) 
    stop("'formula' missing or incorrect")
  if (length(attr(terms(formula[-2L]), "term.labels")) != 1L) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  y <- behrens.test(x = DATA[[1L]], y = DATA[[2L]], ...)
  if (length(y$estimate) == 2L) {
    names(y$estimate) <- paste("mean in group", levels(g))
    names(y$null.value) <- paste("difference in means between",
                                 paste("group", levels(g), collapse = " and "))
  }
  y$data.name <- DNAME
  y
}
