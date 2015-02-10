\name{quantileKM}
\alias{quantileKM}
\title{Compute arbitrary quantile of a Kaplan-Meier estimate of a survival function}
\description{In survival analysis, often some quantiles (mainly the median) of an estimated survival function
are of interest. The survfit function in the 'survival' packages of version \eqn{< 2.35} computed median time to event. However, 
the output of that function was such that the median is not accessible by the user (it can only be read
off the output). This function makes the median and any other quantile accessible by the user.}
\usage{quantileKM(time, event, group = NA, quant = 0.5, conf.level = 0.95, 
    conftype = c("log", "log-log", "plain", "none")[2], conflower = 
    c("usual", "peto", "modified")[1])}
\arguments{
  \item{time}{Event times, censored or observed.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{group}{Quantiles can be computed for several survival curves, defined by \code{group}.}
  \item{quant}{Quantile to be computed. Real number in \eqn{[0, 1]}.}
  \item{conf.level}{Significance level for confidence interval for the time to event quantile.}
  \item{conftype}{Type of confidence interval to be computed. For possible choices see above, and for specifications
  regarding the different options confer the help file of the function \code{survfit}.}
  \item{conflower}{Controls modified lower limits to the curve, the upper limit remains unchanged. See \code{survfit} for details.}
}
\value{
\item{n}{Number of observations used.}
\item{events}{Number of events.}
\item{quantile}{Quantile estimate.}
\item{lower.ci}{Lower limit of confidence interval.}
\item{upper.ci}{Upper limit of confidence interval.}
}

\references{Computation of confidence intervals is done according to

Brookmeyer, R. and Crowley, J. (1982). 
A Confidence Interval for the Median Survival Time. 
\emph{Biometrics}, \bold{38}, 29--41.
}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}

\seealso{Partly based on the function \code{survfit}.}
\examples{
## use Acute Myelogenous Leukemia survival data contained in package 'survival'
time <- leukemia[, 1]; status <- leukemia[, 2]; x <- as.factor(leukemia[, 3])
plot(survfit(Surv(time, status) ~ x, conf.type = "none"), mark = "/", col = 1:2)

## median time to event
qKM <- quantileKM(time, status, group = x, quant = 0.5, conf.level = 0.95, 
    conftype = "log")
qKM

## extract results
qKM$quantities

## comparison to standard function (median time to event not accessible by user)
quant.surv <- survfit(Surv(time, status) ~ x, conf.int = 0.95)
quant.surv

## compute arbitrary quantile
qKM2 <- quantileKM(time, status, group = x, quant = 0.25, conf.level = 0.95, 
    conftype = "log")
qKM2
}

\keyword{htest}
\keyword{survival}






