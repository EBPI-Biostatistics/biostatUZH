#' Plot of population and sample relationship
#' 
#' Plot visualizing the connection between a population
#' and a sample drawn from it.
#'  
#' @param language Language of the text of the plot.
#' @return Generates a plot.
#' @author Kaspar Rufibach \cr \email{kaspar.rufibach@@gmail.com}
#' @keywords dplot aplot
#' @examples
#' populationSamplePlot()
#' @export
populationSamplePlot <- function(language = c("german", "english")){
    stopifnot(!is.null(language))
    language <- match.arg(language)
    
    opar <- par(mar = rep(0,  4))
    on.exit(par(opar))
    
    plot(0, 0, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 0.8), axes = FALSE, type = "n")
    A <- 0; B <- 0.7; C <- 0.1; D <- 0.6; a <- 0.4; b <- 0.6; c <- 0.2; d <- 0.3; epsx <- 0.4
    epsy <- 0.4; e <- 0.05; ee <- 0.15
    
    graphics::polygon(c(A, A, B, B), c(C, D, D, C), col = 5, density = -1)
    graphics::polygon(c(a, a, b, b) + epsx, c(c, d, d, c) + epsy, col = 2, density = -1)
    set.seed(22111905)
    x <- stats::runif(10,  0.1,  0.69)
    y <- stats::runif(10,  0.11,  0.59)
    graphics::points(x,  y,  pch = 19,  col = 2)
    for(i in 1:length(x)){
        graphics::lines(c(x[i], a + epsx), c(y[i], c + epsy))
    }
    graphics::lines(c((A + B)/2, (a + b)/2 + epsx), c(d + epsy + e, d + epsy + e))
    graphics::lines(c((A + B)/2, (A + B)/2), c(D, d + epsy + e))
    graphics::lines(c((a + b)/2 + epsx, (a + b)/2 + epsx), c(d + epsy, d + epsy + e))
    
    if (identical(language, "german")){
        graphics::text((A + B)/2, C-e, expression(paste("Populationsmittelwert ", mu)), cex = 0.8)
        graphics::text((a + b)/2 + epsx, c + epsy-e, "Stichprobe:", cex = 0.8)
        graphics::text((a + b)/2 + epsx, c + epsy-2 * e, expression(paste("Mittelwert ", bar(x))), cex = 0.8)
        graphics::text((A + B)/2-ee, D + 2 * e, "Rueckschluss auf", cex = 0.8)
        graphics::text((A + B)/2-1.2 * ee, D + e, "Populationsmittelwert", cex = 0.8)
        graphics::points((A + B)/2, D, pch = 6)
        graphics::points(a + epsx, c + epsy, pch = 6)
    }
    
    if (identical(language, "english")){
        graphics::text((A + B)/2, C-e, expression(paste("population mean ", mu)), cex = 0.8)
        graphics::text((a + b)/2 + epsx, c + epsy-e, "sample:", cex = 0.8)
        graphics::text((a + b)/2 + epsx, c + epsy-2 * e, expression(paste("mean ", bar(x))), cex = 0.8)
        graphics::text((A + B)/2-ee, D + 2 * e, "draw conclusion for ", cex = 0.8)
        graphics::text((A + B)/2-1.2 * ee, D + e, "population mean", cex = 0.8)
        graphics::points((A + B)/2, D, pch = 6)
        graphics::points(a + epsx, c + epsy, pch = 6)
    }
    invisible(NULL)
}
