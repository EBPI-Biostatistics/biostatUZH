#' Fagan-Nomogram
#' 
#' In
#' 
#' 
#' @aliases faganPlot faganLine
#' @param prob.pre.init Pre-test probabilities to be used for the vertical line
#' to the left. May be a number or a vector withe entries in \eqn{(0,1)}.
#' @param text Size of the text displayed in the plot.
#' @param language Choose language.
#' @param tit Title to be added to plot, or \code{NA} if to be left blank.
#' @param prob.pre A single number or a vector of pre-test probabilities to be
#' drawn a line in the Fagan-plot for.
#' @param lik.ratio Likelihood ratios to be used. Must be of same length as
#' \code{prob.pre}.
#' @return Nothing is returned. \code{faganPlot} draws a Fagan-Nomogram and
#' \code{faganLine} adds lines for specific pre-test probabilities and
#' likelihood ratios to the plot.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Fagan, T.J. (1975). Letter: Nomogram for Bayes Theorem. \emph{N.
#' Engl. J. Med}, \bold{293}, 257.
#' @keywords dplot aplot
#' @examples
#' 
#' # empty Fagan plot
#' faganPlot()
#' 
#' # Fagan lines
#' faganLine(prob.pre = 0.75, lik.ratio = c(0.9/0.37, 0.1/0.63))
#' 
#' @export
faganPlot <- function(prob.pre.init = c(.1, .2, .5, 1, 2, 5, 10, 20, 30, 40, 50, 70, 85), text = 0.8, 
    language = c("german", "english")[1], tit = "Fagan - Nomogramm"){

# draw plot
# pretest probabilities & chances
prob.pre <- prob.pre.init / 100
prob.pre <- sort(c(prob.pre, 0.5, 1 - prob.pre))
cpre <- logit(prob.pre)

# log - likelihood ratios
t <- c(1, 2, 5)
d <- rep(c(1000, 100, 10, 1, 0.1, 0.01), times = 1, each = 3)
rep(t, times = 6)
lr <- c(t / d, 1000)
llr <- log(lr)

# plot parameters
tick <- 0.05  # tick width
text <- 0.7   # text size
par(mar = c(1, 1, 1, 1))

plot(0, 0, type ='n', xlab = '', yaxt = 'n', ylab = '', xaxt = 'n', xlim = 1.3 * c(-1, 1), ylim = c(-8, 7), bty = "n")

# plot pre-test probabilities
segments(-1, min(cpre), -1, max(cpre))
segments(-1 - tick, rev(cpre), -1 + tick, rev(cpre))
text(-1.1, rev(cpre), round(100 * prob.pre, 2), adj = 1, cex = text)

# plot log likelihood ratios
llr <- llr / 2
segments(0, min(llr), 0, max(llr))
segments(-tick, llr, tick, llr)
text(-0.1, llr, round(lr, 3), adj = 1, cex = text)

# plot post-test probabilities
segments(1, min(cpre), 1, max(cpre))
segments(1 - tick, cpre, 1 + tick, cpre)
text(1.1, cpre, round(100 * prob.pre, 2), adj = 0, cex = text)

# other texts
if (language == "english"){
text(-1.3, 0, "Pre-test probability (%)", srt = 90)
text(1.3, 0, "Post-test probability (%)", srt = 90)
text(0, 4.2, "Likelihood ratio")
}

if (language == "german"){
text(-1.3, 0, "Pre-test Wahrscheinlichkeit (%)", srt = 90)
text(1.3, 0, "Post-test Wahrscheinlichkeit (%)", srt = 90)
text(0, 4.2, "Likelihood Quotient")
}

if (identical(tit, NA) == FALSE){mtext(tit, 3, 0)}
}
