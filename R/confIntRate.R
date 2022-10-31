#' Wald and Wilson confidence intervals for Poisson rates
#' 
#' Compute a confidence interval for Poisson rates using the Wald and Wilson methods. 
#' 
#' @aliases confIntRate waldRate wilsonRate
#' @rdname confIntRate
#' @param x Number of events.
#' @param t Total observation time.
#' @param conf.level Confidence level for confidence interval.
#' @return \code{confIntRate} returns a data.frame with columns: \code{x}, \code{t}, \code{type},
#' \code{rate}, \code{lower}, \code{upper}.
#' @author Leonhard Held
#' @seealso \code{\link{confIntRateDiff}}, package \pkg{Hmisc}.
#' @references Held, L., Rufibach, K. and Seifert, B. (2013). Medizinische
#' Statistik - Konzepte, Methoden, Anwendungen. Section 8.2.
#' @examples
#' 
#' confIntRateWald(x = 1, t = 2.2)
#' confIntRateWilson(x = 1, t = 2.2)
#' confIntRate(x = 3:1, t = seq(from = 1, to = 2.2, length.out = 3))
#' 
#' @export
confIntRate <- function(x, t, conf.level = 0.95)
{
    ## input is checked in 'confIntRateWald' and 'confIntRateWilson'
    rbind(data.frame(x = x, t = t, type = "Wald",
                     confIntRateWald(x = x, t = t, conf.level = conf.level)),
          data.frame(x = x, t = t, type = "Wilson",
                     confIntRateWilson(x = x, t = t, conf.level = conf.level)))
}


#' @rdname confIntRate
#' @return \code{confIntRateWald}returns a matrix with columns
#' \code{rate}, \code{lower}, and \code{upper}.
#' @export
confIntRateWald <- function(x, t, conf.level = 0.95) {
    stopifnot(is.numeric(x), length(x) > 0,
              is.finite(x), is.wholenumber(x),
              is.numeric(t), length(t) > 0,
              is.finite(t), t > 0,
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    lambda <- x / t
    q <- stats::qnorm(p = (1 + conf.level) / 2)
    ef <- exp(q / sqrt(x))
    cbind("rate" = lambda, "lower" = lambda / ef, "upper" = lambda * ef)
}

#' @export
waldRate <- function(x, t, conf.level = 0.95){
    .Deprecated(new = "confIntRateWald")
    confIntRateWald(x = x, t = t, conf.level = conf.level)
}


#' @rdname confIntRate
#' @return \code{confIntRateWilson} returns a matrix with columns
#' \code{rate}, \code{lower}, and \code{upper}.
#' @export
confIntRateWilson <- function(x, t, conf.level = 0.95) {
    stopifnot(is.numeric(x), length(x) > 0,
              is.finite(x), is.wholenumber(x),
              is.numeric(t), length(t) > 0,
              is.finite(t), t > 0,
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    lambda <- x / t
    q <- stats::qnorm(p = (1 + conf.level) / 2)
    a <- x + q^2 / 2
    b <- q / 2 * sqrt(4 * x + q^2)
    cbind("rate" = lambda, "lower" = (a - b) / t,
          "upper" = (a + b) / t)
}

#' @export
wilsonRate <- function(x, t, conf.level = 0.95){
    .Deprecated(new = "confIntRateWilson")
    confIntRateWilson(x = x, t = t, conf.level = conf.level)
}


#' Wald and Wilson confidence intervals for Poisson rates differences
#' 
#' Compute a confidence interval for Poisson rate differences using the
#' Wald and Wilson methods. 
#' 
#' @param x vector of length 2, number of events in each group.
#' @param t vector of length 2, total observation time in each group.
#' @param conf.level Confidence level for confidence interval.
#' @return A list with the entries:
#' \item{rd}{Estimated rate difference.}
#' \item{CIs}{Data.frame containing confidence intervals for the rate difference.}
#' @author Leonhard Held
#' @seealso \code{\link{confIntRateWilson}}, \code{\link{confIntRateWald}}
#' @references Held, L., Rufibach, K. and Seifert, B. (2013). Medizinische
#' Statistik - Konzepte, Methoden, Anwendungen. Section 8.2.
#' @keywords htest
#' @examples
#' 
#' confIntRateDiff(x = c(30, 50), t = c(100, 120.1))
#' 
#' @export
confIntRateDiff <- function(x, t, conf.level = 0.95){
    stopifnot(is.numeric(x), length(x) == 2,
              is.finite(x), is.wholenumber(x),
              is.numeric(t), length(t) == 2,
              is.finite(t), t > 0,
              is.numeric(conf.level), length(conf.level) == 1,
              is.finite(conf.level),
              0 < conf.level, conf.level < 1)
    
    ci1 <- confIntRateWilson(x = x[1], t = t[1], conf.level = conf.level)
    ci2 <- confIntRateWilson(x = x[2], t = t[2], conf.level = conf.level)

    diff <- matrix(ci1[1] - ci2[1])
    se.diff <- sqrt(x[1] / t[1]^2 + x[2] / t[2]^2)
    z <- stats::qnorm((1 + conf.level) / 2)
    wald.lower <- diff - z * se.diff
    wald.upper <- diff + z * se.diff
    score.lower <- diff - sqrt((ci1[1] - ci1[2])^2 + (ci2[3] - ci2[1])^2)
    score.upper <- diff + sqrt((ci2[1] - ci2[2])^2 + (ci1[3] - ci1[1])^2)

    result <- matrix(ncol = 2, nrow = 2)
    result[,1] <- c(wald.lower, score.lower)
    result[,2] <- c(wald.upper, score.upper)

    out <- data.frame(type = c("Wald", "Wilson"), result)
    names(out) <- c("type", "lower", "upper")

    list("rd" = diff, "CIs" = out)
}


