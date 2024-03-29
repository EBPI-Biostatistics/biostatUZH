% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tableOR.R
\name{tableOR}
\alias{tableOR}
\title{Results table with odds ratios from logistic regression models for binary or
ordinal data}
\usage{
tableOR(
  model,
  caption = "",
  label = "",
  size = "scriptsize",
  factorNames = NULL,
  table.placement = "ht",
  refLevels = NULL,
  language = c("english", "german"),
  short = FALSE,
  latex = TRUE,
  rmStat = FALSE,
  wald = FALSE
)
}
\arguments{
\item{model}{Object of class \code{"polr"} or \code{"glm"} with family binomial.}

\item{caption}{Character vector of length 1 containing the table caption.}

\item{label}{Character vector of length 1 containing the reference label of
the Latex table.}

\item{size}{Character vector of length 1 setting LaTeX font site, e.g. "small",
"scriptsize", etc.}

\item{factorNames}{A character vector of size k number of factors or
regressor with custom factor labels.}

\item{table.placement}{LaTeX table positioning. The default is "ht".
See \code{?print.xtable} for more information.}

\item{refLevels}{A character vector of size k number of regressors with
custom reference level names. This argument is usuful for objects of class
\code{glm()}, which does not store the reference level labels internally.}

\item{language}{Language of the table, either "english" (default) or "german".}

\item{short}{If \code{TRUE}, factor names are removed from factor
levels. Default is \code{FALSE}.}

\item{latex}{If \code{TRUE} (default) LateX output is produced.}

\item{rmStat}{Logical, if \code{FALSE} (default) output table inludes test
statistics.}

\item{wald}{Logical, if \code{FALSE} (default) Wilson confidence intervals
are computed.}
}
\value{
The table as data.frame and, depending on the value of the argument \code{latex},
the function a print of the LaTeX code.
}
\description{
Result table with odds ratios, 95\%-CI, test statistics, and p-values from
logistic regression models for binary or ordinal variables.
}
\examples{

data <- carData::TitanicSurvival
# relevel: baseline is survived yes.
data$survived <- relevel(x = data$survived, ref = "yes") 
model <- glm(survived ~ sex + age + passengerClass, data = data, family = binomial())
tableOR(model = model, latex = FALSE, short = TRUE,
        refLevels = c("female", "1st"), 
        caption = "Changes in odds for risk of death in the Titanic tragedy.")

## using log regression for ordinal data
data$passengerClass <- factor(x = data$passengerClass, ordered = TRUE)
model <- MASS::polr(passengerClass ~ sex + age, data = data, Hess = TRUE)
tableOR(model = model, latex = FALSE, short = TRUE, 
        caption = "Changes in odds for being in a lower class, i.e. 2nd or 3rd class")

}
\seealso{
\code{\link{tableRegression}}, \code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
}
\author{
Simon Schwab
}
