\name{behrensTest}
\alias{behrensTest}
\title{Computes Behrens version of the t-test for unequal variances.}
\description{
Computes behrens version of the t-test for unequal variances based on approximate
 solution as in Box & Tiao, 1973, Section 2.5.3.
}

\usage{
behrensTest(x, y, conf.level = 0.95, alternative = "two.sided")
}

\arguments{
  \item{x}{A (non-empty) numeric vector of data values.}
  \item{y}{A (non-empty) numeric vector of data values.}
  \item{conf.level}{Confidence level for confidence intervals.}
  \item{alternative}{A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". String must match exactly.}
}

\value{
	\item{lower}{Lower limit of the confidence interval of the difference.}
	\item{lower}{Estimated difference in means.}
	\item{upper}{Upper limit of the confidence interval of the difference.}
	\item{t-value}{t-statistic.}
	\item{df}{Degrees of freedom.}
	\item{p-value}{p-value}
}

\author{.}

\note{
	Based on approximate solution as in Box & Tiao, 1973, Section 2.5.3.
}

\seealso{Armitage, Berry, Matthews, 2002, Section 4.3}



