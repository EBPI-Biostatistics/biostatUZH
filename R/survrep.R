survrep <- function(ftime, fstatus, fgroup = NA, 
        stmt.pl = c("bottomleft", "bottomright", "topright", "subtitle", "none")[1], 
        legend.pl = c("bottomleft", "bottomright", "topright", "none")[3],
        output = c("plain", "text")[2],
        maintitle = "",
        ylbl = "Survival",
        xlbl = "Time", 
        lbls = levels(fgroup),
        dig = 2, conf.level = 0.95){

if (identical(fgroup, NA) == TRUE) {
    fm <- Surv(ftime, fstatus) ~ 1
    np <- 1
    
} else {
    fm <- Surv(ftime, fstatus) ~ fgroup
    np <- length(table(fgroup))
}

out <- match.arg(output, c("plain", "text"))
stmt.placement <- match.arg(stmt.pl, c("bottomleft", "bottomright", "topright", "subtitle",  "none"))
legend.placement <- match.arg(legend.pl, c("bottomleft", "bottomright", "topright", "none"))
if((stmt.placement == legend.placement) & stmt.placement!="none") warning("Legend and Hazard Ratio Statement will be overplotted! Change legend or statment placement option.")


sv <- survfit(fm)
tab <- summary(sv)$table

if(np==1){
medians <- quantileKM(ftime, fstatus, fgroup, conf.level = conf.level)$quantities[3:5]
med.stmt <- paste(round(medians[1], dig), " (", round(medians[2], dig), " - ", round(medians[3], dig), ")", sep = "")
} else {
medians <- quantileKM(ftime, fstatus, fgroup, conf.level = conf.level)$quantities[,3:5]
med.stmt <- paste(round(medians[,1], dig), " (", round(medians[,2], dig), " - ", round(medians[,3], dig), ")", sep = "")
}

plot(sv, col = 1:np, lty = 1:np, conf.int = FALSE, 
        main = maintitle, xlab = xlbl, ylab = ylbl)
summary(sv)
if(np>1){
    cph <- coxph(fm)
    cph.var <- cph$var
    if(np>2) cph.var <- diag(cph.var)
    z <- coef(cph)/sqrt(cph.var)
    if(np==2){
            pval <- 1 - pchisq(survdiff(fm)$chisq, df = length(survdiff(fm)$n) - 1)
    } else  pval <- pnorm(-abs(z))*2
    contra <- paste(lbls[1], " vs. ", lbls[2:np], sep = "")
    hr.txt <- cbind(exp(coef(cph)), exp(confint(cph, level = conf.level)), pval)
    rownames(hr.txt) <- contra

    hr.stmt <- paste(contra, ": HR ", round(hr.txt[,1], dig), " (", round(hr.txt[,2], dig), " - ", round(hr.txt[,3], dig), "), p ", ifelse(hr.txt[,4]<0.0001, "", "= "), format.pval(hr.txt[,4], digits = dig, eps = 0.0001), sep = "") 
    hr.stmt2 <- hr.stmt
    if(!stmt.placement %in% c("topright")) hr.stmt2 <- rev(hr.stmt)
    if(!stmt.placement %in% c("subtitle", "none")){
        mult <- switch(stmt.placement, bottomleft = 1/50, bottomright = 1 - 1/50, topright = 1 - 1/50)
        algn <- switch(stmt.placement, bottomleft = 4, bottomright = 2, topright = 2)
        incr <- (1:(np - 1))/20
        vert <- switch(stmt.placement, bottomleft = 0 + incr, bottomright = 0 + incr, topright = 1 - incr, none = 0)
        text(x = max(sv$time)*mult, y = vert, labels = hr.stmt2, pos = algn)
    }
    if(stmt.placement=="subtitle"){    
        mtext(hr.stmt2, side = 3, line = 0:(np - 2))
    }
    if(legend.placement!="none") legend(legend.placement, lbls, col = 1:np, lty = 1:np, bty = "n")
}

if(out == "plain"){
    tmp.out <- list(med = medians)
    if(np>1) tmp.out$hr <- hr.txt
} else {
    tmp.out <- list(med = med.stmt)
    if(np>1) tmp.out$hr <- hr.stmt
}

return(tmp.out)

}
