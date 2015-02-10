\name{tableRegression}
\alias{tableRegression}
\title{Present regression output}
\description{Modifies regression output in a desired way.}
\usage{
tableRegression(model, stats = NULL,
                col.nam = NULL, row.nam = NULL, intercept = NULL,
                text = "english", text.ci = text, 
                eps.pvalue = 0.0001, digits = NULL, big.mark = "'",
                xtable = TRUE, align = NULL, caption = NULL, label = NULL)
}
\arguments{
  \item{model}{Object of class \code{lm},  \code{glm}, \code{coxph}, or
    \code{list} (for Weibull models).}
  \item{stats}{character vector with stats, choose from \code{"estimate",
      "exp.estimate", "standarderror", "t.value", "ci.95", "p.value"}.}
  \item{col.nam}{character vector, same length and order as in
    \code{stats}. If the percentage sign is used, two
    backslashes are needed.}
  \item{row.nam}{Names of rows, character vector}
  \item{intercept}{logical, should intercept be provided or not. If
    model is a binomial regression, intercept is set \code{FALSE}. If intercept is set \code{TRUE}, the first line of the
    summary output is removed (intercepts are not available for weibull or cox models,
    because they do not provide any intercept value).}
  \item{text}{Character, german or english, used e.g. for naming the
    \var{intercept}.}
  \item{text.ci}{Character, german or english or none, referring to the language used for the
    confidence interval, see \code{\link[reporttools]{displayCI}}.}
  \item{eps.pvalue}{p-values smaller than \code{eps.pvalue} will be
    formatted as '< eps.pvalue'.}
  \item{digits}{Vector of length \code{stats}, digits used for each column.}
  
  \item{big.mark}{Character vector, see \code{\link{format}}.}
  \item{align}{Character vector, see \code{\link[xtable]{xtable}}}
  \item{xtable}{Logical, presented as \code{xtable} or just matrix?}
  \item{caption}{Character vector, see \code{\link[xtable]{xtable}}}
  \item{label}{Character vector, see \code{\link[xtable]{xtable}}}
}

\value{Depending on the value of the argument \code{xtable} a
  matrix or a \code{xtable} object is returned.}
\details{In \code{stats}: 
  \itemize{
    \item If \code{t.value} is choosen, the
    \code{z.value} might be taken, depending on the model.
    \item For lm-models: \code{ci.95} calculates a confidence interval
    for the estimate.
    \item For glm- anx coxph-models: \code{ci.95} calculates a confidence interval
    for the exp(estimate).
  }}

\author{Sina Rueeger\cr \email{sina.rueeger@uzh.ch}}
\seealso{\code{\link[xtable]{xtable}}}
\note{If the percentage sign is used, two
    backslashes are needed (not only one as in latex). Because of
    problems when running the package, the two backslashes and the
    percent sign are omitted
    in the examples.}

\examples{

## Linear model
## ---------------
mod.lm <- lm(Sepal.Length ~ Sepal.Width, data = iris)
mod.lm1 <- lm(Sepal.Length ~ .^2, data = iris) 

## using default style
tableRegression(mod.lm)

## choosing columns, plus columns and row naming german
tableRegression(mod.lm1, stats = c("estimate", "t.value", "p.value"), text = "german")

## adapt row names, plus special format for ci
tableRegression(mod.lm, row.nam = c("Intercept", "Width Sepal"), text.ci = "none")

## Poission model
## (example from ?glm)
## --------------
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
mod.glm.pois <- glm(counts ~ outcome + treatment, family=poisson())
tableRegression(mod.glm.pois)


## Logistic model
## -------------
dat <- survival::rats
dat$rx <- factor(dat$rx, labels = c(" (A)", " (B)"))
mod.glm.bin <- glm(status ~ litter + rx, family = binomial, data
= dat)

tableRegression(mod.glm.bin,
                stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
                text = "english", digits = rep(3, 5),
                caption = "Here goes the caption.", label = "mod:logit")

## including intercept
tableRegression(mod.glm.bin,
                stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
                text = "english", digits = rep(3, 5),
                caption = "Here goes the caption.", label = "mod:logit",
                intercept = TRUE)


## Cox model
## (example from ?survival::coxph)
## -------------
dat <- list(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 

library("survival")
mod.cox <- coxph(Surv(time, status) ~ x, dat)
mod.cox1 <- coxph(Surv(time, status) ~ x + factor(sex), dat)
mod.cox2 <- coxph(Surv(time, status) ~ x + strata(sex), dat)

tableRegression(mod.cox)
tableRegression(mod.cox1)
tableRegression(mod.cox2)


## Weibull
## (example from biostatUZH::WeibullReg)
## -------------
data(larynx)
mod.wb <- WeibullReg(Surv(time, death) ~ factor(stage) + age, data=larynx)
tableRegression(mod.wb)

}
