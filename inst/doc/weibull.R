### R code from vignette source 'weibull.Rnw'

###################################################
### code chunk number 1: source
###################################################
library(biostatUZH)
data(larynx)
larynx$stage2 <- as.numeric(larynx$stage == 2)
larynx$stage3 <- as.numeric(larynx$stage == 3)
larynx$stage4 <- as.numeric(larynx$stage == 4)


###################################################
### code chunk number 2: WeibullReg
###################################################
WeibullReg(Surv(time, death) ~ stage2 + stage3 + stage4 + age, data=larynx)


###################################################
### code chunk number 3: coxph
###################################################
summary(coxph(Surv(time, death) ~ stage2 + stage3 + stage4 + age, data=larynx))$conf.int


###################################################
### code chunk number 4: DiagPlot
###################################################
fm <- Surv(time, death) ~ stage
diagWR <- WeibullDiag(fm, larynx, labels=c("Stage I", "Stage II", "Stage III", "Stage IV"))


