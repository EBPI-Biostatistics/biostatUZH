
## computes behrens version of the t-test for unequal variances
## based on approximate solution as in Box & Tiao, 1973, Section 2.5.3
## see also Armitage, Berry, Matthews, 2002, Section 4.3

behrensTest <- function(x, y, conf.level = 0.95, alternative="two.sided"){
    alpha <- 1-conf.level
    m1 <- mean(x)
    s1 <- sd(x)
    n1 <- sum(!is.na(x))
    m2 <- mean(y)
    s2 <- sd(y)
    n2 <- sum(!is.na(y))
  
    nu1 <- n1-1
    nu2 <- n2-1
    diff <- m1-m2
    diff.se <- sqrt(s1^2/n1+s2^2/n2)
    
    cos2.phi <- (s2^2/n2)/(diff.se^2)
    sin2.phi <- (s1^2/n1)/(diff.se^2)
    
    f1  <- nu2/(nu2-2)*cos2.phi + nu1/(nu1-2)*sin2.phi 
    f2 <- nu2^2/((nu2-2)^2*(nu2-4))*cos2.phi^2 + nu1^2/((nu1-2)^2*(nu1-4))*sin2.phi^2 
    df <- 4+f1^2/f2
    a <- sqrt((df-2)/df*f1)
    t.value <- diff/(a*diff.se)
    if(alternative=="two.sided"){
        c <- qt(alpha/2, df=df, lower.tail = FALSE)
        lower <- diff-a*diff.se*c
        upper <- diff+a*diff.se*c
        p <-  2*pt(abs(t.value), df=df, lower.tail = FALSE)
    }
    if(alternative=="less"){
        c <- qt(alpha, df=df, lower.tail = FALSE)
        lower <- -Inf
        upper <- diff+a*diff.se*c
        p <-  pt(t.value, df=df, lower.tail = TRUE)
    }
    if(alternative=="greater"){
        c <- qt(alpha, df=df, lower.tail = FALSE)
        lower <- diff-a*diff.se*c
        upper <- Inf
        p <-  pt(t.value, df=df, lower.tail = FALSE)
    }
    res <- c("lower"=lower, "diff"=diff, "upper"=upper, "t-value"=t.value, "df"=df, "p-value"=p)
    return(res)
}
