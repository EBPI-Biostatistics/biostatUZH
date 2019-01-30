\name{tableOR}
\alias{tableOR}
\title{Results table with odds ratios from logistic regression models for binary or ordinal data}
\description{Result table with odds ratios, 95\%-CI, test statistics, and p-values from logistic regression models for binary or ordinal variables.}

\usage{
tableOR(model, caption="", label="", size="scriptsize", factorNames=NULL,
                   table.placement = "ht", refLevels=NULL, lang="english", short=FALSE, latex=TRUE, rmStat=FALSE)
}

\arguments{
  \item{model}{an object of class \code{"polr"} or \code{"glm"}, the latter has to be from the family binomial.}
  \item{caption}{Table caption.}
  \item{label}{A string containing the LaTeX table reference label.}
  \item{size}{A string to set LaTeX font site, e.g. small, scriptsize, etc.}
  \item{factorNames}{A character vector of size k number of factors or regressor with custom factor labels.}
  \item{table.placement}{LaTeX table positioning.}
    \item{refLevels}{A character vector of size k number of regressors with custom reference level names. This option is designed when using \code{glm()} which does not store the reference level labels internally.}
  \item{lang}{Language of the confidence intervals term, "english"" (default) or "german".}
  \item{short}{A logical, if \code{TRUE} factor names are removed from factor levels. Default is \code{FALSE}.}
  \item{latex}{Logical, if \code{TRUE} (default) LateX output is produced.}
  \item{rmStat}{Logical, if \code{FALSE} (default) output table inludes test statistics.}
}

\value{
  Depending on the value of the argument \code{latex}, the function
  either prints LaTeX code or returns a data frame.
}

\author{Simon Schwab}

\seealso{\code{\link[xtable]{xtable}}}

\examples{
dat = carData::TitanicSurvival
dat$survived = relevel(dat$survived, ref = "yes") # relevel: baseline is survived yes.
model = glm(survived ~ sex + age + passengerClass,
            data = dat, family = binomial())
labels = c("female", "1st") # reference levels of the two categorial variables
tableOR(model, latex = FALSE, short = TRUE, refLevels = labels, caption = "Changes in odds for risk of death in the Titanic tragedy.")

## using log regression for ordinal data
dat$passengerClass = factor(dat$passengerClass, ordered = TRUE)
model = MASS::polr(passengerClass ~ sex + age, data = dat, Hess = TRUE)
tableOR(model, latex = FALSE, short = TRUE, caption = "Changes in odds for being in a lower class, i.e. 2nd or 3rd class")
}
