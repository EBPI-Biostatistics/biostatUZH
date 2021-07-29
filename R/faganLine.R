#' @export
faganLine <- function(prob.pre, lik.ratio = c(1, 1)){

# draw line: test is positiv
logits.post <- log(lik.ratio[1]) + logit(prob.pre)
prob.post <- exp(logits.post) / (1 + exp(logits.post))
segments(-1, logit(1 - prob.pre), 1, logits.post, lwd = 2, col = 2)

# draw line: test is negativ
logits.post <- log(lik.ratio[2]) + logit(prob.pre)
prob.post <- exp(logits.post) / (1 + exp(logits.post))
segments(-1, logit(1 - prob.pre), 1, logits.post, lwd = 2, col = 3)
}



