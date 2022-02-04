\name{confIntFisherTrafo}
\Rdversion{1.1}
\alias{confIntFisherTrafo}
\title{Confidence interval for correlation coefficient using Fisher's transformation}
\description{Compute a confidence interval for a correlation coefficient \eqn{r} using the variance-stabilizing
transformation

\deqn{z = \tanh^{-1}(r) = 0.5 \log((1 + r) / (1 - r)),}

known as Fisher's \eqn{z}-transformation. By means of this transformation, \eqn{r} is approximately normally distributed
with variance \eqn{(n-3)^{-1}} independent of the true correlation \eqn{\rho}, enabling construction of a Wald-type confidence interval. Back-transformation yields
a confidence interval for the correlation coefficient. An advantage of this approach is that the limits of 
the confidence interval are contained in \eqn{(-1, 1)}.
}
\usage{confIntFisherTrafo(var1, var2, pp = c(0.025, 0.975), meth = 
    "spearman", type = "t")}
\arguments{
  \item{var1}{Vector containing first variable.}
  \item{var2}{Vector containing first variable.}  
  \item{pp}{Vector in \eqn{R^2} that contains \eqn{\alpha / 2} and \eqn{1 - \alpha/2}, where
  \eqn{alpha} is the confidence level of the confidence interval.}
  \item{meth}{Correlation coefficient to be used: \code{pearson} or \code{spearman}.}
  \item{type}{Quantile to be used: \code{z} or \code{t}.}
}
\value{Yields a list with entries:
  \item{estimate}{Value of correlation coefficient.}
  \item{ci}{Computed confidence interval.}
  \item{p.value}{\eqn{p}-value for a test on \eqn{\rho = 0} based on the transformation.}
  \item{n}{Number of observations.}
  \item{p2}{\eqn{p}-value based on the \code{R} function \code{cor.est}.}
}
\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}
\examples{
n <- 40
x <- runif(n)
y <- 2 * x + 0.5 * rnorm(n)
plot(x, y)
confIntFisherTrafo(x, y, pp = c(0.025, 0.975), meth = "spearman", type = "t")
}
\keyword{htest}
