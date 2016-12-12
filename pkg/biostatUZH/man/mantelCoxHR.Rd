\name{mantelCoxHR}
\alias{mantelCoxHR}
\title{Mantel-Cox estimator of the hazard ratio}
\description{Estimate the hazard ratio and compute a confidence interval for it via a special application of the 
Mantel-Haenszel procedure. A separate 2 x 2 table is constructed for each event time. The underlying assumption is 
that the hazard ratio is constant over the follow-up period.}
\usage{mantelCoxHR(time, event, group, conf.level = 0.95)}
\arguments{
  \item{time}{Event times, censored or observed.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{group}{Factor with two levels, e.g. treatment group.}
  \item{conf.level}{Significance level for confidence interval for hazard ratio.}
}

\value{
\item{mantelCox.hr}{Hazard ratio estimate.}
\item{ci.hr}{Wald confidence interval at the level specified by \code{conf.level}.}
\item{p.val.logrank}{\eqn{p}-value of logrank test for a comparison of the survival curves between the two groups.}
\item{Cox-reg}{Hazard ratio estimated via Cox-regression.}
}

\references{
Kirkwood, B.R. and Sterne, J.A.C. (2003). \emph{Essential Medical Statistics}. Blackwell Science. See p. 283 ff.
}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com} }

\note{Note that in general the Mantel-Cox estimate and the hazard ratio estimate received via Cox-regression
\emph{do not coincide}.}

\examples{
## use Acute Myelogenous Leukemia survival data contained in package 'survival'
time <- leukemia[, 1]; status <- leukemia[, 2]; x <- as.factor(leukemia[, 3])
mantelCoxHR(time, status, x)
}

\keyword{htest}
\keyword{survival}
