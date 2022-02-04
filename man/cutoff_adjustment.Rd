\name{cutoff_adjustment}
\alias{cutoff_adjustment}
\title{Calculated adjusted p-value and adjusted hazard ratio when using minimal p-value approach.}
\description{
Adjustement of the p-value is done according to Hilsenbeck & Clark 1996 "Practical p-value adjustment for optimally selected cutpoints".
Adjustment of the hazard ratio is done via the heuristic shrinkage proposed by van Houwelingen and Le Cessie "Predictive value of statistical models", 1990, Statistics in medicine / Schumacher et al. 2012, pp. 424f..
}

\usage{
cutoff_adjustment(time, time2, event, outcome, trim = 2,
                  na.omit = FALSE, data, firth = FALSE, reverse = FALSE)
}

\arguments{
  \item{time}{A character string, corresponds to 'time' from survival::Surv}
  \item{time2}{A character string, corresponds to 'time2' from survival::Surv}
  \item{event}{A character string, corresponds to 'event' from survival::Surv}
  \item{outcome}{A character string for the continuous outcome variable.}
  \item{trim}{Number of events that should be excluded on the lower and upper end 
  to avoid complete separation.}
  \item{na.omit}{Logical value. If TRUE, rows with missing values will be excluded.
  Default is FALSE.}
  \item{data}{Data frame.}
  \item{firth}{Logical value. If TRUE, then Firth Cox-Regression will be used. 
  Default is FALSE. Package "coxphf" must be installed.}
  \item{reverse}{Logical value. If TRUE, then the event variable will be switched within 
  the cox regression. Default is FALSE.}
  
}

\value{
	\item{n.event}{Number of events.}
	\item{meas.event}{Values of variable outcome associated with events.}
	\item{bound.l}{Lower bound of truncation, i.e. smallest value considered after trim.}
	\item{bound.u}{upper bound of truncation, i.e. largest value considered after trim.}
	\item{meas.event.trunc}{Truncated set of values of covariates.}
	\item{cutpoints}{Set of cutpoints used in the minimal p-value approach.}
	\item{p.unadj.min}{Minimal p-value.}
	\item{cutoff}{Estimated cut-off value.}
	\item{p.adj}{Adjusted p-value.}
	\item{c.hat}{Estimated shrinkage coefficient.}
	\item{beta.adj}{Adjusted beta coefficient. Direction depends on argument "reverse".}
	\item{ci.beta.adj}{95\% confidence interval for the adjusted beta coefficient. Direction depends on argument "reverse".}
	\item{hr.adj}{Adjusted hazard ratio. Direction depends on argument "reverse".}
	\item{ci.hr.adj}{95\% confidence interval for the adjusted hazard ratio. Direction depends on argument "reverse".}
}

\author{Niels Hagenbuch, Klaus Steigmiller}

\seealso{Houwelingen and Le Cessie "Predictive value of statistical models", 1990, Statistics in medicine / Schumacher et al. 2012, pp. 424f in "Handbook of statistics in clinical oncology" / Hilsenbeck & Clark 1996 "Practical p-value adjustment for optimally selected cutpoints"}

\examples{
\dontrun{
cutoff_adjustment(time = "time.num", time2 = "end.time", event = "event", 
                  na.omit = TRUE,
                  outcome = "lactat", trim = 2,
                  data = ddsurv, 
                  firth = FALSE, reverse = TRUE)
}
}


