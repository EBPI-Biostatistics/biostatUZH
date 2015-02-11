
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
