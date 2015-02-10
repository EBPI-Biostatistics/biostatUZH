\name{survrep}
\alias{survrep}
\title{Survival Analysis Results for Reports}
\description{
  The analysis of survival data most often requires the
  production of Kaplan-Meier plots, estimated median survival time, and
  univariate estimates of hazard ratios (with associated confidence
  intervals and p-values). Often, this is performed separately for
  several possible risk factors. This function simplifies these tasks by
  providing a unified interface to simultaneously produce graphics,
  estimate median survival time, and compute hazard ratios with
  associated statistics from the Cox proportional hazards model.
}

\usage{
survrep(ftime, fstatus, fgroup = NA,
        stmt.pl = c("bottomleft", "bottomright", "topright", "subtitle", "none")[1],
        legend.pl = c("bottomleft", "bottomright", "topright", "none")[3],
        output = c("plain", "text")[2],
        maintitle = "",
        ylbl = "Survival",
        xlbl = "Time",
        lbls = levels(fgroup),
        dig = 2, 
        conf.level = 0.95)
}
\arguments{
  \item{ftime}{Event times, censored or observed.}
  \item{fstatus}{Censoring indicator, 1 for event, 0 for censored.}
  \item{fgroup}{Kaplan-Meier plots, median survival and hazard rates
    will be computed by \code{fgroup}. If \code{fgroup} is missing,
    Kaplan-Meier plots and median survival will be produced for the
    entire dataset.}
  \item{stmt.pl}{Placement of the statement giving hazard ratios,
    confidence intervals and p-values in the plot. Ignored if
    \code{fgroup = NA}.}
  \item{legend.pl}{Placement of the legend in the plot. Ignored if
    \code{fgroup = NA}. A warning will be produced if \code{stmt.pl =
      legend.pl}.}
  \item{output}{Type of output, ``plain'' gives plaintext output for
    median survival and hazard ratios suitable for the creation of
    tables etc., while ``text''  gives output suitable for the text of
    reports. ``text''-type output is always printed in the plots unless
    \code{stmt.pl="none"}.}
  \item{maintitle}{Main title of the plot. Corresponds to the
    \code{main} option of \code{plot.survfit}.}
  \item{ylbl}{Label for the y-axis of the plot. Corresponds to the
    \code{ylab} option of \code{plot.survfit}.}
  \item{xlbl}{Label for the x-axis of the plot. Corresponds to the
    \code{xlab} option of \code{plot.survfit}.}
  \item{lbls}{Labels used to produce the legend, as well as the
    contrasts printed in the hazard ratio statements.}
  \item{dig}{Number of significant digits used in rounding for
    ``text''-type output.}
  \item{conf.level}{Significance level for confidence interval for the
    median survival and hazard ratio estimates.}
}

\value{
  The Kaplan-Meier plot is printed, as well as a list containing:
  \item{med}{Matrix of median survival, and corresponding confidence
    intervals. If \code{output = "text"}, this will be a vector of text
    giving median survival estimates and confidence intervals.}
  \item{events}{Matrix of hazard ratios, corresponding confidence
    intervals, and p-values from the Wald test. If \code{output = "text"},
    this will be a vector of text.}
}

\author{Sarah R. Haile}

\seealso{\code{survfit}, \code{coxph}, \code{quantileKM}}

\examples{
## use Acute Myelogenous Leukemia survival data contained in package 'survival'
time <- leukemia[, 1]; status <- leukemia[, 2]; x <- leukemia[, 3]
survrep(time, status, NA)
survrep(time, status, x)
survrep(time, status, x, stmt.pl = "bottomright", out = "plain", main = "Acute Myeloid Leukemia")
survrep(time, status, x, stmt.pl = "bottomright", out = "text", main = "AML")
survrep(time, status, x, stmt.pl = "subtitle", out = "text", main = "Acute Myeloid Leukemia")

## using Larynx data contained in package 'survival'
## to show behavior with risk factors with >2 levels
data(larynx)
larynx$stage <- factor(larynx$stage, 1:4,  c("I", "II", "III", "IV"))
with(larynx, survrep(time, death, stage, out = "t",
                     stmt.pl = "bottomright", main = "Larynx"))
with(larynx, survrep(time, death, stage, out = "t",
                     stmt.pl = "subtitle", main = ""))
with(larynx, survrep(time, death, stage, out = "p",
                     stmt.pl = "topright", legend.pl = "bottomright"))
}

\keyword{survival}
