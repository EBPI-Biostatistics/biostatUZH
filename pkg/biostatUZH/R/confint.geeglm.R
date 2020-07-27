
# # From package broom version 0.5.3 (not exported from broom Namespace)
# in order to avoid external dependences
# 
# broom:::confint.geeglm

confint.geeglm.broom <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-mult*Std.err,
                      upr=Estimate+mult*Std.err))
  rownames(citab) <- rownames(cc)
  citab[parm,]
}