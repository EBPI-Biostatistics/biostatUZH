
## computes confidence interval for the difference in AUC of two paired tests
## cases is a m x 2 matrix
## controls is a n x 2 matrix
## each row corresponds to measurements from one case/control

confIntIndependentAUCDiff <- function(casesA, controlsA, casesB, controlsB,
                                      type="Wald", conf.level = 0.95)
{
    ncontrolsA <- nrow(controlsA)
    ncasesA <- nrow(casesA)
    ncontrolsB <- nrow(controlsB)
    ncasesB <- nrow(casesB)
    auc <- numeric()

    resA <- confIntAUC(casesA, controlsA, conf.level = conf.level)
    resB <- confIntAUC(casesB, controlsB, conf.level = conf.level)
    factor <- qnorm((1 + conf.level) / 2)

    type <- match.arg(type, c("Wald", "Logit"))

    if(type=="Wald"){
        ## take intervals on original scale
        aucA <- resA[1,3]
        aucB <- resB[1,3]
        D <- aucA - aucB
        lowerA <- resA[1,2]
        lowerB <- resB[1,2]
        upperA <- resA[1,4]
        upperB <- resB[1,4]
        ## compute and combine standard errors
        seA <- (upperA-lowerA)/(2*factor)
        seB <- (upperB-lowerB)/(2*factor)
        se.D <- sqrt(seA^2+seB^2)
        D.lower <- D - factor * se.D
        D.upper <- D + factor * se.D

        res <- data.frame(matrix(NA, ncol = 4))
        colnames(res) <- c("outcome", "lower", "estimate", "upper")
        res[1, 2:4] <- resA[1,2:4] # Wald interval on original scale
        res[2, 2:4] <- resB[1,2:4] # Wald interval on original scale
        res[3, 2:4] <- c(D.lower, D, D.upper)
        res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    }

    if(type=="Logit"){
        ## take intervals on logit scale
        aucA <- resA[2,3]
        aucB <- resB[2,3]
        D <- aucA - aucB
        lowerA <- resA[2,2]
        lowerB <- resB[2,2]
        upperA <- resA[2,4]
        upperB <- resB[2,4]

        ## Apply Newcombe trick
        D.lower <- D - sqrt((aucA - lowerA) ^ 2 + (aucB - upperB) ^ 2)
        D.upper <- D + sqrt((aucA - upperA) ^ 2 + (aucB - lowerB) ^ 2)

        res <- data.frame(matrix(NA, ncol = 4))
        colnames(res) <- c("outcome", "lower", "estimate", "upper")
        res[1, 2:4] <- resA[2,2:4] # Wald interval on logit scale
        res[2, 2:4] <- resB[2,2:4] # Wald interval on logit scale
        res[3, 2:4] <- c(D.lower, D, D.upper)
        res[, 1] <- c("AUC Test 1", "AUC Test 2", "AUC Difference")
    }

    return(res)
}
