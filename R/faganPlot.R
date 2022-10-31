#' Fagan-Nomogram
#'
#' Draw the Fagan-Nomogram.
#' 
#' @rdname faganPlot
#' @aliases faganPlot faganLine
#' @param probPreInit Pre-test probabilities in percent to be used for the vertical line
#' to the left. A number or a vector with entries in \eqn{(0,100)}.
#' @param cex Text size.
#' @param language Either "english" (default) or "german".
#' @param title Title of the plot.
#' @return \code{faganPlot} draws a Fagan-Nomogram. 
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @references Fagan, T.J. (1975). Letter: Nomogram for Bayes Theorem. \emph{N.
#' Engl. J. Med}, \bold{293}, 257.
#' @keywords dplot aplot
#' @examples
#' 
#' # empty Fagan plot
#' faganPlot()
#' 
#' # add Fagan lines
#' faganLine(probPre = 0.75, likRatio = c(0.9 / 0.37, 0.1 / 0.63))
#' 
#' @export
faganPlot <- function(probPreInit = c(.1, .2, .5, 1, 2, 5, 10, 20, 30, 40, 50, 70, 85),
                      cex = 0.7, language = c("english", "german"),
                      title = "Fagan - Nomogramm"){
    
    stopifnot(is.numeric(probPreInit), length(probPreInit) > 0,
              is.finite(probPreInit), 0 < probPreInit, probPreInit < 100,
              is.numeric(cex), length(cex) == 1,
              is.finite(cex), cex > 0,
              !is.null(language))
    language <- match.arg(language)
    stopifnot(is.character(title), length(title) == 1)
    
# draw plot
# pretest probabilities & chances
    probPre <- probPreInit / 100
    probPre <- sort(c(probPre, 0.5, 1 - probPre))
    cpre <- stats::qlogis(probPre)
    
# log - likelihood ratios
    t <- c(1, 2, 5)
    d <- rep(c(1000, 100, 10, 1, 0.1, 0.01), times = 1, each = 3)
    #rep(t, times = 6) # Commented this out because it does nothing
    lr <- c(t / d, 1000)
    llr <- log(lr)
    
# plot parameters
    tick <- 0.05  # tick width
    oldPar <- graphics::par(mar = c(1, 1, 1, 1))
    on.exit(oldPar)
    plot(0, 0, type ='n', xlab = '', yaxt = 'n', ylab = '', xaxt = 'n',
         xlim = 1.3 * c(-1, 1), ylim = c(-8, 7), bty = "n")
    
                                        # plot pre-test probabilities
    graphics::segments(-1, min(cpre), -1, max(cpre))
    graphics::segments(-1 - tick, rev(cpre), -1 + tick, rev(cpre))
    graphics::text(-1.1, rev(cpre), round(100 * probPre, 2), adj = 1, cex = cex)
    
                                        # plot log likelihood ratios
    llr <- llr / 2
    graphics::segments(0, min(llr), 0, max(llr))
    graphics::segments(-tick, llr, tick, llr)
    graphics::text(-0.1, llr, round(lr, 3), adj = 1, cex = cex)
    
                                        # plot post-test probabilities
    graphics::segments(1, min(cpre), 1, max(cpre))
    graphics::segments(1 - tick, cpre, 1 + tick, cpre)
    graphics::text(1.1, cpre, round(100 * probPre, 2), adj = 0, cex = cex)
    
                                        # other texts
    if (language == "english"){
      graphics::text(-1.3, 0, "Pre-test probability (%)", srt = 90)
      graphics::text(1.3, 0, "Post-test probability (%)", srt = 90)
      graphics::text(0, 4.2, "Likelihood ratio")
    }    
    if (language == "german"){
      graphics::text(-1.3, 0, "Pre-test Wahrscheinlichkeit (%)", srt = 90)
      graphics::text(1.3, 0, "Post-test Wahrscheinlichkeit (%)", srt = 90)
      graphics::text(0, 4.2, "Likelihood Quotient")
    }
    
    if (identical(title, NA) == FALSE){
      graphics::mtext(title, 3, 0)
    }
    invisible(NULL)
}



#' @rdname faganPlot 
#' @param probPre Numric vector specifying the pre-test probability.
#' @param likRatio Numeric vector specifying the likelihood ratio.
#' @return \code{faganLine} adds lines to a plot obtained by \code{faganPlot}.
#' @export
faganLine <- function(probPre, likRatio = c(1, 1)){

                                        # draw line: test is positiv
    logitsPost <- log(likRatio[1]) + stats::qlogis(probPre)
    probPost <- exp(logitsPost) / (1 + exp(logitsPost))
    graphics::segments(-1, stats::qlogis(1 - probPre), 1, logitsPost, lwd = 2, col = 2)
    
                                        # draw line: test is negativ
    logitsPost <- log(likRatio[2]) + stats::qlogis(probPre)
    probPost <- exp(logitsPost) / (1 + exp(logitsPost))
    graphics::segments(-1, stats::qlogis(1 - probPre), 1, logitsPost, lwd = 2, col = 3)
    invisible(NULL)
}
