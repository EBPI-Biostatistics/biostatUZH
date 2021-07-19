## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function NumEvents


#' Compute number of events for a survival endpoint
#' 
#' Calculates either the required total number of events or the power.
#' 
#' 
#' @param HR Hazard ratio of experimental group vs. control group.
#' @param sig.level Significance level.
#' @param power Desired power.
#' @param n.events Required total number of events.
#' @param alloc.ratio Allocation ratio: Ratio of the number of patients in the
#' experimental group divided by the number of patients in the control group.
#' @param non.inf.margin Non-inferiority margin.
#' @param type Study type. Either "sup" or "superiority" for superiority
#' studies or "noninf" or "non-inferiority" for non-inferiority studies.
#' Default is "sup".
#' @param alternative In c("one.sided","two.sided"), depending on whether the
#' alternative hypothesis is one- or two-sided.
#' @return Returns either the required total number of events or the power.
#' @author Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@@uzh.ch}
#' @references Collett, D. (2015). \emph{Modelling Survival Data in Medical
#' Research}.  New York: Chapman and Hall/CRC.
#' 
#' Schoenfeld, D.A. (1983). Sample-Size Formula for the Proportional-Hazards
#' Regression Model.  \emph{Biometrics}, \bold{39}, 499--503.
#' @examples
#' 
#' NumEvents(HR = 0.65, sig.level = 0.05, power = 0.9,
#'           alloc.ratio = 1, type = "sup", alternative = "two.sided")
#' 
NumEvents <- function(HR, sig.level=0.05, power=NULL, n.events=NULL,
                     alloc.ratio=1, non.inf.margin=NULL, type="sup",
                     alternative="two.sided"){
    stopifnot(HR>0)
    logHR <- log(HR)
    stopifnot(alternative %in% c("two.sided","one.sided"))
    tails <- ifelse(alternative=="one.sided",1,2)
    z_alpha <- qnorm(sig.level/tails)
    p <- alloc.ratio/(alloc.ratio+1)
    p.factor <- p*(1-p)
    ## calculate # events
    if(is.null(n.events)){
        z_beta <- qnorm(1-power)
        c <- (z_alpha+z_beta)^2
        ## superiority case
        if(type %in% c("superiority","sup")){
            d <- c/(p.factor*logHR^2)
            return(ceiling(d))
        }
        ## non-inferiority case
        if(type %in% c("non-inferiority","noninf")){
            stopifnot(!is.null(non.inf.margin))
            d <- c/(p.factor*(logHR - non.inf.margin)^2)
            return(ceiling(d))
        }
    }
    ## calculate power
    if(!is.null(n.events)){
        ## superiority case
        if(type %in% c("superiority","sup")){
            power <- 1 - pnorm(logHR*sqrt(n.events*p.factor)-z_alpha)
            return(power)
        }
        ## non-inferiority case
        if(type %in% c("non-inferiority","noninf")){
            stopifnot(!is.null(non.inf.margin))
            power <- 1 - pnorm((logHR-non.inf.margin)*sqrt(n.events*p.factor)-z_alpha)
            return(power)
        }
    }
    stop("Wrong distribution")
}
