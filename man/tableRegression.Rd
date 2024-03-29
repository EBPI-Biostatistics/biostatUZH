% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tableRegression.R
\name{tableRegression}
\alias{tableRegression}
\title{Format coefficient tables of regression models}
\usage{
tableRegression(
  model,
  stats = NULL,
  col.names = NULL,
  row.names = NULL,
  intercept = NULL,
  text = c("english", "german"),
  text.ci = text,
  eps.pvalue = 1e-04,
  digits = NULL,
  strict = FALSE,
  big.mark = "'",
  xtable = TRUE,
  align = NULL,
  caption = NULL,
  label = NULL,
  vars = NULL,
  ...
)
}
\arguments{
\item{model}{Object of class \code{\link[stats]{lm}},
\code{\link[stats]{glm}}, \code{negbin} (obtained by
\code{\link[MASS]{glm.nb}}), \code{coxph} (obtained by
\code{\link[survival]{coxph}}), and list obtained by \code{\link{Weibull}}.}

\item{stats}{character vector with stats chosen from "estimate",
"exp.estimate", "standarderror", "t.value", "ci.95", and "p.value".}

\item{col.names}{Character vector of same length and order as in \code{stats}.
A percentage sign must be escaped by two backslashes.}

\item{row.names}{Character vector of row names.}

\item{intercept}{Logical vector of length one indicating whether to provide
an intercept or not. If intercept is set TRUE, the first line of the summary
output is removed. If the model is a binomial regression, intercept is set
FALSE. Intercepts are not available for Weibull or Cox models, because they
do not provide any intercept value.}

\item{text}{Either "english" (default) or "german" indicating the used
language names.}

\item{text.ci}{Either "english", "german" or "none". The language used to
denote confidence interval, see \code{\link{formatCI}}.}

\item{eps.pvalue}{If \code{strict = FALSE}, p-values smaller than 
\code{eps.pvalue} will be formatted as "< eps.pvalue". Otherwise this
argument is ignored.}

\item{digits}{Vector of length \code{stats}, digits used for each column.}

\item{strict}{Either \code{TRUE} or \code{FALSE} (default). If \code{FALSE}, 
p-values are formatted with \code{\link[biostatUZH]{formatPval}}. Otherwise
\code{\link[biostatUZH]{formatPvalStrict}} is used to format p-values. 
This argument is thus only relevant when the argument \code{stats} is 
either \code{NULL} or contains \code{"p.value"}.}

\item{big.mark}{Character vector as in \code{\link[base]{format}}.}

\item{xtable}{If TRUE, a Latex table is returned, otherwise a data.frame is
returned.}

\item{align}{See \code{\link[xtable]{xtable}}.}

\item{caption}{See \code{\link[xtable]{xtable}}.}

\item{label}{See \code{\link[xtable]{xtable}}.}

\item{vars}{Specify the variables for which regression summaries should be
printed.  The argument \code{vars} takes a string vector with the names of
the coefficients in the model.}

\item{...}{Arguments passed to \code{\link{print.xtable}}.}
}
\value{
Depending on the value of the \code{xtable} argument, the function
either prints and returns LaTeX code representing the produced table of
coefficients, or it returns the corresponding data frame.
}
\description{
Formats output from the regression model functions: \code{\link[stats]{lm}},
\code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}},
\code{\link[survival]{coxph}}, and \code{\link{Weibull}}.
}
\details{
In \code{stats}: \itemize{ \item If \code{t.value} is chosen, the
\code{z.value} might be taken, depending on the model.  \item For lm-models:
\code{ci.95} calculates a confidence interval for the estimate.  \item For
glm- and coxph-models: \code{ci.95} calculates a confidence interval for the
exp(estimate). }
}
\examples{

## Linear model
## ---------------
mod.lm <- lm(Sepal.Length ~ Sepal.Width, data = iris)
mod.lm1 <- lm(Sepal.Length ~ .^2, data = iris) 

tableRegression(model = mod.lm)

## strict argument is used to force p-values to have a certain number
## of digits (here 4).
tableRegression(model = mod.lm, digits = c(2, 2, 4), strict = TRUE)
tableRegression(model = mod.lm, digits = c(2, 2, 4), strict = FALSE)

## choosing columns, columns and row naming in german
tableRegression(model = mod.lm1, stats = c("estimate", "t.value", "p.value"),
                text = "german")

## adapt row names, plus special format for ci
tableRegression(model = mod.lm, row.names = c("Intercept", "Width Sepal"),
                text.ci = "none")

## Poisson model
## (example from ?glm)
## --------------
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
mod.glm.pois <- glm(counts ~ outcome + treatment, family=poisson())
tableRegression(model = mod.glm.pois)


## Negative binomial model
## --------------
mod.glm.nb <- glm.nb(Days ~ Sex + Age, data = quine)
tableRegression(model = mod.glm.nb,
    caption = paste("NegBin model. Estimated dispersion:",
        sprintf("\%4.2f ($se=\%4.2f$).", mod.glm.nb$theta, mod.glm.nb$SE.theta)),
    label = "tab:glm.nb")



## Logistic model
## -------------
dat <- survival::rats
dat$rx <- factor(dat$rx, labels = c(" (A)", " (B)"))
mod.glm.bin <- glm(status ~ litter + rx, family = binomial, data = dat)

tableRegression(model = mod.glm.bin,
                stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
                text = "english", digits = rep(3, 5),
                caption = "Here goes the caption.", label = "mod:logit")

## including intercept
tableRegression(model = mod.glm.bin,
                stats = c("estimate", "exp.estimate", "ci.95", "t.value", "p.value"),
                text = "english", digits = rep(3, 5),
                caption = "Here goes the caption.", label = "mod:logit",
                intercept = TRUE)


## Cox model
## (example from ?survival::coxph)
## -------------
dat <- list(time = c(4, 3, 1, 1, 2, 2, 3), 
            status = c(1, 1, 1, 0, 1, 1, 0), 
            x = c(0, 2, 1, 1, 1, 0, 0), 
            sex = c(0, 0, 0, 0, 1, 1, 1)) 

mod.cox <- coxph(Surv(time, status) ~ x, data = dat)
mod.cox1 <- coxph(Surv(time, status) ~ x + factor(sex), data = dat)
mod.cox2 <- coxph(Surv(time, status) ~ x + strata(sex), data = dat)

tableRegression(model = mod.cox)
tableRegression(model = mod.cox1)
tableRegression(model = mod.cox2)


## Weibull
## (example from WeibullReg)
## -------------
data("larynx")
mod.wb <- weibullReg(Surv(time, death) ~ factor(stage) + age,
                     data = larynx)
tableRegression(model = mod.wb)

}
\seealso{
\code{\link{tableOR}}, \code{\link[xtable]{xtable}},
\code{\link[stats]{lm}}, \code{\link[stats]{glm}},
\code{\link[MASS]{glm.nb}} \code{\link[survival]{coxph}},
\code{\link{Weibull}}.
}
\author{
Sina Rueeger with contributions by Sebastian Meyer.
}
