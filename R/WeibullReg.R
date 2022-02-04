WeibullReg <-
function (formula, data = parent.frame(), conf.level = 0.95) 
{
    m <- survreg(formula, data, dist = "weibull")
    mle <- ConvertWeibull(m, conf.level)
    return(list(formula = formula, coef = mle$vars, HR = mle$HR, 
        ETR = mle$ETR, summary = summary(m)))
}
