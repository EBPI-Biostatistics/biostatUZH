#' Confidence interval for correlation coefficient using Fisher's
#' transformation
#' 
#' Compute a confidence interval for a correlation coefficient \eqn{r} using
#' the variance-stabilizing transformation
#' 
#' \deqn{z = \tanh^{-1}(r) = 0.5 \log((1 + r) / (1 - r)),}
#' 
#' known as Fisher's \eqn{z}-transformation. By means of this transformation,
#' \eqn{r} is approximately normally distributed with variance \eqn{(n-3)^{-1}}
#' independent of the true correlation \eqn{\rho}, enabling construction of a
#' Wald-type confidence interval. Back-transformation yields a confidence
#' interval for the correlation coefficient. An advantage of this approach is
#' that the limits of the confidence interval are contained in \eqn{(-1, 1)}.
#' 
#' 
#' @param var1 Vector containing first variable.
#' @param var2 Vector containing first variable.
#' @param pp Vector in \eqn{R^2} that contains \eqn{\alpha / 2} and \eqn{1 -
#' \alpha/2}, where \eqn{alpha} is the confidence level of the confidence
#' interval.
#' @param meth Correlation coefficient to be used: \code{pearson} or
#' \code{spearman}.
#' @param type Quantile to be used: \code{z} or \code{t}.
#' @return Yields a list with entries: \item{estimate}{Value of correlation
#' coefficient.} \item{ci}{Computed confidence interval.}
#' \item{p.value}{\eqn{p}-value for a test on \eqn{\rho = 0} based on the
#' transformation.} \item{n}{Number of observations.} \item{p2}{\eqn{p}-value
#' based on the \code{R} function \code{cor.est}.}
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @keywords htest
#' @examples
#' 
#' n <- 40
#' x <- runif(n)
#' y <- 2 * x + 0.5 * rnorm(n)
#' plot(x, y)
#' confIntFisherTrafo(x, y, pp = c(0.025, 0.975), meth = "spearman", type = "t")
#'
#' @export
confIntFisherTrafo <- function(var1, var2, pp = c(0.025, 0.975), meth = "spearman", type = "t"){

# default values
estimate <- NA
ci1 <- NA
p.val <- NA
p2 <- NA

x <- cbind(var1, var2)
x <- x[is.na(apply(x, 1, sum)) == 0, ]
x <- matrix(x, ncol = 2)
n <- nrow(x)

if (n >= 2){
    corEst <- cor.test(x[, 1], x[, 2], use = "pairwise", method = meth)
    estimate <- corEst$estimate
    p2 <- corEst$p.value
    }

if (n >= 4){
    #cor.est <- cor(x, use = "pairwise", method = meth)[1, 2]
    rho <- 0.5 * log((1 + corEst$estimate) / (1 - corEst$estimate))

    if (type == "z"){
        quant <- qnorm(pp)
        p.val <- pnorm(abs(rho) * sqrt(n - 3))}
    if (type == "t"){
        quant <- qt(pp, df = n - 3)
        p.val <- pt(abs(rho) * sqrt(n - 3), df = n - 3)}

    ci1 <- rho + quant / sqrt(n - 3)
    ci1 <- (exp(2 * ci1) - 1) / (exp(2 * ci1) + 1)

    # does the same
    # ci1 <- tanh(atanh(corEst$estimate) + quant / sqrt(n - 3))

    p.val <- 2 * (1 - p.val)
    }

res <- list("estimate" = as.numeric(estimate), "ci" = ci1, "p.value" = as.numeric(p.val), "n" = n, "p2" = p2)
return(res)
}








