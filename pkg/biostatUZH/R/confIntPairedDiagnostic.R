
# assumes data for Diseased and nonDiseased are given in a 2x2 table each
# with negative diagostic test counts in the first column and first row
# and positive diagostic test counts in the second column and first row
# This is as in Pepe (2003), Section 3.3

confIntPairedDiagnostic <- function(Diseased, nonDiseased, conf.level = 0.95, adjust=FALSE)
{
    stopifnot(is.wholenumber(Diseased), is.wholenumber(nonDiseased),
              Diseased>0, nonDiseased>0,  conf.level<1,
              conf.level>0)
    stopifnot(nrow(Diseased)==2, ncol(Diseased)==2, nrow(nonDiseased)==2, ncol(nonDiseased)==2)

    ratio <- function(x) return(x[1]/x[2])
    rPF <- colSums(Diseased)/rowSums(Diseased)
    rNF <- colSums(nonDiseased)/rowSums(nonDiseased)
    rLR <- rPF/rNF
    
    se.log.rTPF <- sqrt((sum(Diseased)-sum(diag(Diseased)))/(rowSums(Diseased)[2]*colSums(Diseased)[2]))
    se.log.rFPF <- sqrt((sum(nonDiseased)-sum(diag(nonDiseased)))/(rowSums(nonDiseased)[2]*colSums(nonDiseased)[2]))
    se.log.rFNF <- sqrt((sum(Diseased)-sum(diag(Diseased)))/(rowSums(Diseased)[1]*colSums(Diseased)[1]))
    se.log.rTNF <- sqrt((sum(nonDiseased)-sum(diag(nonDiseased)))/(rowSums(nonDiseased)[1]*colSums(nonDiseased)[1]))
    se.log.rLRplus <- sqrt(se.log.rTPF^2 + se.log.rFPF^2)
    se.log.rLRminus <- sqrt(se.log.rFNF^2 + se.log.rTNF^2)
    
    rEstimates <- c(rPF[2], rNF[1], rPF[1], rNF[2], rev(rLR))
    if(adjust)
        conf.level <- sqrt(conf.level)
    z <- qnorm((1 + conf.level) / 2)
    EF <- exp(z*c(se.log.rTPF, se.log.rTNF, se.log.rFNF, se.log.rFPF, se.log.rLRplus, se.log.rLRminus))

    res <- data.frame(matrix(NA, ncol=4, nrow=6))
    colnames(res) <- c("type", "lower", "estimate", "upper")
    res[, 1] <- c("rTPF", "rTNF", "rFNF", "rFPF", "rLRplus", "rLRminus")
    res[, 2] <- rEstimates/EF
    res[, 3] <- rEstimates
    res[, 4] <- rEstimates*EF
    
    return(res)
}
