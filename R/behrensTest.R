#' Behrens version of the t-test for unequal variances
#'
#' Behrens version of the t-test for unequal variances.
#' @param x Numeric vector of data values.
#' @param y Numeric vector of data values.
#' @param alternative Character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater", and "less".
#' @param conf.level Confidence level of the interval, default is 0.95.
#' @details Computes behrens version of the t-test for unequal variances
#' based on approximate solution as in Box & Tiao, 1973, Section 2.5.3;
#' see also Armitage, Berry, Matthews, 2002, Section 4.3.
#' @examples
#' weightWomen <- c(38.9, 61.2, 73.3, 21.8, 63.4, 64.6, 48.4, 48.8, 48.5)
#' weightMen <-   c(67.8, 60  , 63.4, 76,   89.4, 73.3, 67.3, 61.3)
#'
#' behrensTest(x = weightWomen, y = weightMen)
#' behrensTest(x = weightWomen, y = weightMen, alternative="less")
#' @export
behrensTest <- function(x, y, conf.level = 0.95,
                        alternative = c("two.sided", "greater", "less")){
    stopifnot(is.numeric(x), length(x) >= 1, is.finite(x),
              is.numeric(y), length(y) >= 1, is.finite(y),
              is.numeric(conf.level), length(conf.level) == 1,
              0 < conf.level, conf.level < 1,
              !is.null(alternative))
    alternative <- match.arg(alternative)
    
    alpha <- 1 - conf.level
    m1 <- mean(x)
    s1 <- sd(x)
    n1 <- sum(!is.na(x))
    m2 <- mean(y)
    s2 <- sd(y)
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
    t.value <- diff / (a * diff.se)
    
    if(alternative == "two.sided"){
        c <- qt(alpha / 2, df = df, lower.tail = FALSE)
        lower <- diff - a * diff.se * c
        upper <- diff + a * diff.se * c
        p <- 2 * pt(abs(t.value), df = df, lower.tail = FALSE)
    }
    if(alternative == "less"){
        c <- qt(alpha, df = df, lower.tail = FALSE)
        lower <- -Inf
        upper <- diff + a * diff.se * c
        p <-  pt(t.value, df = df, lower.tail = TRUE)
    }
    if(alternative == "greater"){
        c <- qt(alpha, df = df, lower.tail = FALSE)
        lower <- diff - a * diff.se * c
        upper <- Inf
        p <-  pt(t.value, df = df, lower.tail = FALSE)
    }
    c("lower" = lower, "diff" = diff, "upper" = upper,
      "t-value" = t.value, "df" = df, "p-value" = p)
}
