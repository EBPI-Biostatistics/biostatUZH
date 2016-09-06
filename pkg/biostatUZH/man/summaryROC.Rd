\name{summaryROC}
\alias{summaryROC}
\title{Function to compute ROC curve and an asymptotic confidence interval for AUC}
\description{This function computes for values of a continuous variable of a group of cases and a group of controls
the ROC curve. Additionally, the AUC including an asymptotic confidence interval is provided.}
\usage{summaryROC(cases, controls, conf.level = 0.95)}
\arguments{
  \item{cases}{Values of the continuous variable for the cases.}
  \item{controls}{Values of the continuous variable for the controls.}
  \item{conf.level}{Confidence level of confidence interval.}
}
\value{
  \item{x.val}{1-specificity of the test, so the values on the \eqn{x}-axis of a ROC plot.}
  \item{y.val}{Sensitivity of the test, so the values on the \eqn{y}-axis of a ROC plot.}
  \item{ppvs}{Positive predictive values for each cutoff.}
  \item{npvs}{Negative predictive values for each cutoff.}
  \item{cutoffs}{Cutoffs used (basically the pooled marker values of cases and controls).}
  \item{res.mat}{Collects the above quantities in a matrix, including Wilson confidence intervals, computed at at confidence level \code{conf.level}.}
  \item{auc}{Area under the ROC curve. This is equal to the value of the Mann-Whitney test statistic.}
  \item{auc.var}{Variance of AUC.}
  \item{auc.var.norm}{Variance of AUC if data is assumed to come from a bivariate normal distribution.}
  \item{lowCI}{Lower limit of Wald confidence interval.}
  \item{upCI}{Upper limit of Wald confidence interval.}
  \item{logitLowCI}{Lower limit of a Wald confidence interval received on logit scale.}
  \item{logitUpCI}{Upper limit of a Wald confidence interval received on logit scale.}
}

\references{The original reference for computation of confidence intervals is:

Hanley, J.A. and McNeil, B.J. (1982).
The meaning and use of the area under the curve.
\emph{Radiology}, \bold{143}, 29--36.

See also:

Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press.
}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}.
Part of the function was derived from code by Andrea Riebler.}

\note{The confidence intervals are only valid if observations are \emph{independent}.}

\seealso{Similar functionality is provided in the package \pkg{ROCR}. However, this latter package offers no computation
of confidence intervals.}

\examples{

set.seed(1977)
ns <- c(50, 40)
truth <- c(rep(0, ns[1]), rep(1, ns[2]))
estimates <- c(rnorm(ns[1]), rnorm(ns[2], mean = 0.5, sd = 1.5))
cases <- estimates[truth == 1]
controls <- estimates[truth == 0]
res <- summaryROC(cases, controls, conf.level = 0.95)

# display results
res
res$res.mat

# plot ROC curve
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
    xlab = "1 - specificity", ylab = "sensitivity", pty = 's')
segments(0, 0, 1, 1, lty = 2)
lines(res$x.val, res$y.val, type = 'l', col = 2, lwd = 2, lty = 2)
}

\keyword{htest}

