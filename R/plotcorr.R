################################################################################
### Part of the R package "biostatUZH".
### Free software under the terms of the GNU General Public License (version 2
### or later) a copy of which is available at http://www.R-project.org/Licenses
###
### A fork of ellipse::plotcorr() by Duncan Murdoch
###
### Copyright (C) 2015 Sebastian Meyer
################################################################################



#' Compact Visualization of a Correlation Matrix
#' 
#' This function plots a correlation matrix using ellipse-shaped glyphs for
#' each entry.  The ellipse represents a level curve of the density of a
#' bivariate normal with the matching correlation.
#' 
#' This is a fork of the original \code{\link[ellipse]{plotcorr}} function from
#' the \pkg{ellipse} package as at version 0.3-8. The arguments \code{numbers},
#' \code{type}, and \code{diag} have been replaced by \code{lower.panel},
#' \code{upper.panel}, and \code{diag.panel} similar to \code{\link{pairs}}.
#' This enables displaying numbers in one triangle and ellipses in the other.
#' However, there is no support for \code{diag = FALSE} in the original sense
#' of the function.
#' 
#' The ellipses being plotted will be tangent to a unit character square, with
#' the shape chosen to match the required correlation.
#' 
#' @param corr A matrix containing entries between \code{-1} and \code{1} to be
#' plotted as correlations.
#' @param outline Whether the ellipses should be outlined in the default
#' colour.
#' @param col Which colour(s) to use to fill the ellipses (recycled to
#' \code{length(corr)}). The special setting \code{col = TRUE} will fill the
#' ellipses according to the color range \code{colorRampPalette(c("blue",
#' "white", "red"))(11)} mapped to the correlations as shown in the example
#' below.
#' @param lower.panel,upper.panel,diag.panel each panel can be either
#' \code{NULL}, \code{"ellipse"}, or \code{"number"}. \code{NULL} disables the
#' respective part, and \code{"number"} plots numerical correlations in place
#' of ellipses. Correlations will be rounded to a single decimal place. The
#' \pkg{ellipse} package is required for \code{"ellipse"} panels.
#' @param bty,axes,xlab,ylab,asp,mar,cex.lab,... Graphical parameters which
#' will be passed to \code{\link{plot}} when plotting.
#' @param cex Graphical parameter which will be passed to \code{\link{text}}
#' when plotting.
#' @author (of this fork) Sebastian Meyer, (of the original version) Duncan
#' Murdoch
#' @seealso \code{\link[ellipse]{ellipse}}
#' @references Murdoch, D.J. and Chow, E.D. (1996). A graphical display of
#' large correlation matrices. The American Statistician 50, 178-180.
#' @keywords hplot
#' @examples
#' 
#' if (requireNamespace("ellipse")) {
#'   ## Plot the correlation matrix for the mtcars data full model fit 
#'   data("mtcars")
#'   fit <- lm(mpg ~ ., mtcars)
#'   corr.fit <- summary(fit, correlation = TRUE)$correlation
#'   plotcorr(corr.fit, col = "gray")
#' 
#'   ## with default color coding
#'   plotcorr(corr.fit, col = TRUE)
#' 
#'   ## Colour the ellipses and order by correlations with miles/gallon
#'   corr.mtcars <- cor(mtcars)
#'   ord <- order(corr.mtcars[1,])
#'   xc <- corr.mtcars[ord, ord]
#'   colors <- colorRampPalette(c("blue", "white", "red"))(11)
#'   plotcorr(xc, col = colors[5*xc + 6])
#' }
#' 
"plotcorr" <-
  function (corr, outline = TRUE, col = TRUE,
            lower.panel = "ellipse", upper.panel = "number", diag.panel = NULL,
            bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
            cex.lab = par("cex.lab"), cex = 0.75*par("cex"), mar = 0.1 + c(2,2,4,2), ...)
{
    savepar <- par(pty = "s", mar = mar)
    on.exit(par(savepar))

    if (is.null(corr)) return(invisible())
    if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) 
			   || (round(max(corr, na.rm = TRUE), 6) > 1)) 
	stop("Need a correlation matrix")

    plot.new()
    par(new = TRUE)

    rowdim <- dim(corr)[1]
    coldim <- dim(corr)[2]

    rowlabs <- dimnames(corr)[[1]]
    collabs <- dimnames(corr)[[2]]
    if (is.null(rowlabs)) rowlabs <- 1:rowdim
    if (is.null(collabs)) collabs <- 1:coldim
    rowlabs <- as.character(rowlabs)
    collabs <- as.character(collabs)

    col <- if (isTRUE(col)) {
        colorRampPalette(c("blue", "white", "red"))(11)[5*corr + 6]
    } else {
        rep(col, length.out = length(corr))
    }
    dim(col) <- dim(corr)

    if (!is.null(lower.panel)) {
        lower.panel <- match.arg(lower.panel, choices = c("ellipse", "number"))
    }
    if (!is.null(upper.panel)) {
        upper.panel <- match.arg(upper.panel, choices = c("ellipse", "number"))
    }
    if (!is.null(diag.panel)) {
        diag.panel <- match.arg(diag.panel, choices = c("ellipse", "number"))
    }
    
    cols <- 1:coldim
    rows <- 1:rowdim
    
    maxdim <- max(length(rows), length(cols))

    plt <- par('plt')
    xlabwidth <- max(strwidth(rowlabs[rows],units='figure',cex=cex.lab))/(plt[2]-plt[1])
    xlabwidth <- xlabwidth*maxdim/(1-xlabwidth)
    ylabwidth <- max(strwidth(collabs[cols],units='figure',cex=cex.lab))/(plt[4]-plt[3])
    ylabwidth <- ylabwidth*maxdim/(1-ylabwidth)

    plot(c(-xlabwidth-0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), 
	 type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, 
	 cex.lab = cex.lab, ...)
    text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
    text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], 
	 srt = 90, adj = 0, cex = cex.lab)
    mtext(xlab,1,0)
    mtext(ylab,2,0) 
    mat <- diag(c(1, 1))
    plotcorrInternal <- function()
    {
      panel <- if (i == j) diag.panel else if (i > j) lower.panel else upper.panel
      if (is.null(panel)) return()
      if (panel == "ellipse") {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse::ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline) lines(ell)
      } else {
        text(j + 0.3, length(rows) + 1 - i, round(10 * corr[i, j], 0),
             adj = 1, cex = cex)
      }
    }
    for (i in 1:dim(corr)[1]) {
      for (j in 1:dim(corr)[2]) {
          plotcorrInternal()
      }
    }
    invisible()
}
