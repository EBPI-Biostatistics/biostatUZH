printResults <- function(theta, se.theta, conf.level=0.95, digits=3){
  z <- qnorm((1+conf.level)/2)
  ci.l <- theta - z*se.theta
  ci.u <- theta + z*se.theta
  p <-  2*pnorm(abs(theta)/se.theta, lower.tail=FALSE)
  effect <- round(theta, digits=digits)
  ci <- formatCI(c(ci.l, ci.u), digits=digits, text="english")
  p <- formatPval(p)
  res <- c(effect, ci, p)
  names(res) <- list("Effect", 
                     sprintf('%d%% Confidence Interval', conf.level*100),
                     "P-value")
  print(t(res), quote=FALSE)
}
