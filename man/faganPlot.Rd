\name{faganPlot}
\alias{faganPlot}
\alias{faganLine}
\title{Fagan-Nomogram}
\description{In }
\usage{
faganPlot(prob.pre.init = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 
    20, 30, 40, 50, 70, 85), 
    text = 0.8, language = c("german", "english")[1], 
    tit = "Fagan - Nomogramm")
faganLine(prob.pre, lik.ratio = c(1, 1))}
\arguments{
  \item{prob.pre.init}{Pre-test probabilities to be used for the vertical line to the left. May be a number or a vector withe entries in \eqn{(0,1)}.}
  \item{text}{Size of the text displayed in the plot.}
  \item{language}{Choose language.}
  \item{tit}{Title to be added to plot, or \code{NA} if to be left blank.}
  \item{prob.pre}{A single number or a vector of pre-test probabilities to be drawn a line in the Fagan-plot for.}
  \item{lik.ratio}{Likelihood ratios to be used. Must be of same length as \code{prob.pre}.}
}
\value{Nothing is returned. \code{faganPlot} draws a Fagan-Nomogram and \code{faganLine} adds lines for specific
pre-test probabilities and likelihood ratios to the plot.}

\references{
Fagan, T.J. (1975).
Letter: Nomogram for Bayes Theorem.
\emph{N. Engl. J. Med}, \bold{293}, 257.
}

\author{Kaspar Rufibach \cr \email{kaspar.rufibach@gmail.com}}

\examples{
# empty Fagan plot
faganPlot()

# Fagan lines
faganLine(prob.pre = 0.75, lik.ratio = c(0.9/0.37, 0.1/0.63))
}

\keyword{dplot}
\keyword{aplot}
