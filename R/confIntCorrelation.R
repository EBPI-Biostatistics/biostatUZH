#' Confidence interval for correlation coefficient using Fisher's
#' transformation
#' 
#' Computes a confidence interval for a correlation coefficient \eqn{r} using
#' the variance-stabilizing transformation
#' \deqn{z = \tanh^{-1}(r) = 0.5 \log((1 + r) / (1 - r)),}
#' known as Fisher's \eqn{z}-transformation.
#' Independent of the true correlation \eqn{\rho}, \eqn{z} is approximately
#' normally distributed with variance \eqn{(n-3)^{-1}}. This enables the construction
#' of a Wald-type confidence interval. Back-transformating this interval yields a confidence
#' interval for \eqn{r}. An advantage of this method is
#' that the interval is contained in \eqn{(-1, 1)}.
#' 
#' @param x Vector containing the first variable.
#' @param y Vector of same length as \code{x} containing the second variable.
#' @param conf.level Confidence level for confidence interval. Default is 0.95.
#' @param method Correlation coefficient to be used: "spearman" (default) or "pearson".
#' @param type Quantile to be used: "t" (default) or "z".
#' @return List with entries:
#' \item{estimate}{Value of correlation coefficient.}
#' \item{ci}{Computed confidence interval.}
#' \item{p.value}{\eqn{p}-value for a test on \eqn{\rho = 0} based on the transformation.}
#' \item{n}{Number of observations.}
#' \item{p2}{\eqn{p}-value based on the \code{R} function \code{cor.est}.}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @keywords htest
#' @examples
#' 
#' n <- 40
#' x <- runif(n = n)
#' y <- 2 * x + 0.5 * rnorm(n = n)
#' confIntCorrelation(x = x, y = y)
#' 
#' @export
confIntCorrelation <- function(x, y, conf.level = 0.95,
                               method = c("spearman", "pearson"),
                               type = c("t", "z")){

    stopifnot(is.numeric(x), length(x) >= 2L, is.finite(x),
              is.numeric(y), length(y) == length(x), is.finite(y),
              length(conf.level) == 1L, is.finite(conf.level),
              0 < conf.level, conf.level < 1,
              !is.null(method))
    method <- match.arg(method)
    stopifnot(!is.null(type))
    type <- match.arg(type)

    ## default values
    ci <- NA_real_
    p.value <- NA_real_

    n <- length(x)
    alpha <- 1 - conf.level
    corEst <- stats::cor.test(x = x, y = y, use = "pairwise", method = method,
                              conf.level = conf.level)
    estimate <- corEst$estimate
    p2 <- corEst$p.value
    
    if (n >= 4L){
        rho <- 0.5 * log((1 + corEst$estimate) / (1 - corEst$estimate))
        
        if (type == "z"){
            quant <- stats::qnorm(1 - alpha / 2)
            p.value <- stats::pnorm(abs(rho) * sqrt(n - 3))
        } else if (type == "t"){
            quant <- stats::qt(1 - alpha / 2, df = n - 3)
            p.value <- stats::pt(abs(rho) * sqrt(n - 3), df = n - 3)
        }
        
        ci <- rho + c(-1, 1) * quant / sqrt(n - 3)
        ci <- (exp(2 * ci) - 1) / (exp(2 * ci) + 1)
        ## does the same
        ## ci <- tanh(atanh(corEst$estimate) + quant / sqrt(n - 3))
        p.value <- 2 * (1 - p.value)
    }

    list("estimate" = as.numeric(estimate), "ci" = ci,
         "p.value" = as.numeric(p.value), "n" = n, "p2" = p2)
}

#' Deprecated: use confIntCorrelation instead
#' @param var1 Vector containing first variable.
#' @param var2 Vector containing first variable.
#' @param pp Vector in \eqn{R^2} that contains \eqn{\alpha / 2} and \eqn{1 - \alpha/2},
#' where \eqn{alpha} is the confidence level of the confidence interval.
#' @param meth Correlation coefficient to be used: \code{pearson} or \code{spearman}.
#' @param type Quantile to be used: \code{z} or \code{t}.
#' @return List with entries:
#' \item{estimate}{Value of correlation coefficient.}
#' \item{ci}{Computed confidence interval.}
#' \item{p.value}{\eqn{p}-value for a test on \eqn{\rho = 0} based on the transformation.}
#' \item{n}{Number of observations.}
#' \item{p2}{\eqn{p}-value based on the \code{R} function \code{cor.est}.}
#' @export
confIntFisherTrafo <- function(var1, var2, pp = c(0.025, 0.975), meth = "spearman", type = "t"){
    .Deprecated("confIntCorrelation")
    confIntCorrelation(x = var1, y = var2, conf.level = pp[2L], method = meth, type = type)
}








