#' Compute number of events or power for a survival endpoint
#'
#' For a given hazard ratio, significance level, and power,
#' the function returns the required total number of events. \cr
#'
#' For a given hazard ratio, significance level, and number of events,
#' the function returns the corresponding power.
#'
#' @param HR Hazard ratio of experimental group vs. control group.
#' @param sig.level Significance level with default 0.05.
#' @param power Desired power.
#' @param n.events Required total number of events.
#' @param alloc.ratio Allocation ratio: Ratio of the number of patients in the
#' experimental group divided by the number of patients in the control group.
#' @param non.inf.margin Non-inferiority margin. Default is 0 and implies a
#' superiority trial.
#' @param alternative "two.sided" (default) or "one.sided". Specifies whether
#' the alternative hypothesis is one- or two-sided.
#' @return If \code{power} is specified, the function returns the required total
#' number of events. If \code{n.events} is specified, the function returns the
#' power.
#' @author Uriah Daugaard and Leonhard Held \cr \email{leonhard.held@@uzh.ch}
#' @references Collett, D. (2015). \emph{Modelling Survival Data in Medical
#' Research}.  New York: Chapman and Hall/CRC.
#'
#' Schoenfeld, D.A. (1983). Sample-Size Formula for the Proportional-Hazards
#' Regression Model. \emph{Biometrics}, \bold{39}, 499--503.
#' @note Replaces \code{\link{NumEvents}}.
#' @examples
#'
#' survEvents(HR = 0.65, sig.level = 0.05, power = 0.9, alloc.ratio = 1,
#'            non.inf.margin = 0, alternative = "two.sided")
#'
#' survEvents(HR = 0.65, sig.level = 0.05, power = 0.9, alloc.ratio = 1,
#'            non.inf.margin = .5, alternative = "one.sided")
#'
#' survEvents(HR = 0.65, sig.level = 0.05, n.events = 35, alloc.ratio = 1,
#'            non.inf.margin = 0, alternative = "two.sided")
#'
#' survEvents(HR = 0.65, sig.level = 0.05, n.events = 35, alloc.ratio = 1,
#'            non.inf.margin = 1, alternative = "one.sided")

#'
#' @export
survEvents <- function(HR, sig.level = 0.05, power = NULL, n.events = NULL,
                       alloc.ratio = 1, non.inf.margin = 0,
                       alternative = c("two.sided", "one.sided")) {

    stopifnot(is.numeric(HR),
              length(HR) == 1,
              is.finite(HR),
              HR > 0,

              is.numeric(sig.level),
              length(sig.level) == 1,
              is.finite(sig.level),
              0 < sig.level, sig.level < 1,

              sum(c(is.null(power), is.null(n.events))) == 1)
    if (!is.null(power)) {
        stopifnot(length(power) == 1,
                  is.finite(power),
                  0 < power, power < 1)
    } else {
        stopifnot(length(n.events) == 1,
                  is.finite(n.events),
                  0 < n.events)
    }
    stopifnot(is.numeric(alloc.ratio),
              length(alloc.ratio) == 1,
              is.finite(alloc.ratio),
              0 < alloc.ratio,

              is.numeric(non.inf.margin),
              length(non.inf.margin) == 1,
              is.finite(non.inf.margin),
              0 <= non.inf.margin,

              !is.null(alternative))
    alternative <- match.arg(alternative)

    logHR <- log(HR)
    z_alpha <- if (alternative == "one.sided") qnorm(sig.level) else qnorm(sig.level / 2)
    p <- alloc.ratio / (alloc.ratio + 1)
    p.factor <- p * (1 - p)
    ## calculate # events
    if (is.null(n.events)) {
        z_beta <- qnorm(1 - power)
        c <- (z_alpha + z_beta)^2
        d <- c / (p.factor * (logHR - non.inf.margin)^2)
        return(ceiling(d))
    } else {
        ## superiority case
        power <- 1 - pnorm((logHR - non.inf.margin) * sqrt(n.events * p.factor) - z_alpha)
        return(power)
    }
    stop("Function should not arrive here")
}


## code originally written by Uriah Daugaard, slightly edited by Leonhard Held

### Function NumEvents


#' Compute number of events for a survival endpoint
#'
#' This funtion is deprecated. Use \code{\link{survEvents}} instead.
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
#' @seealso \code{\link{survEvents}}
#' @examples
#'
#' ## suppress deprecated warning
#' suppressWarnings(NumEvents(HR = 0.65, sig.level = 0.05, power = 0.9,
#'                            alloc.ratio = 1, type = "sup",
#'                            alternative = "two.sided"))
#'
#' @export
NumEvents <- function(HR, sig.level = 0.05, power = NULL, n.events = NULL,
                     alloc.ratio = 1, non.inf.margin = NULL, type = "sup",
                     alternative = "two.sided") {
    .Deprecated(new = "survEvents")
    stopifnot(HR > 0)
    logHR <- log(HR)
    stopifnot(alternative %in% c("two.sided", "one.sided"))
    tails <- ifelse(alternative == "one.sided", 1, 2)
    z_alpha <- qnorm(sig.level / tails)
    p <- alloc.ratio / (alloc.ratio + 1)
    p.factor <- p * (1 - p)
    ## calculate # events
    if (is.null(n.events)) {
        z_beta <- qnorm(1 - power)
        c <- (z_alpha + z_beta)^2
        ## superiority case
        if (type %in% c("superiority", "sup")) {
            d <- c / (p.factor * logHR^2)
            return(ceiling(d))
        }
        ## non-inferiority case
        if (type %in% c("non-inferiority", "noninf")) {
            stopifnot(!is.null(non.inf.margin))
            d <- c / (p.factor * (logHR - non.inf.margin)^2)
            return(ceiling(d))
        }
    }
    ## calculate power
    if (!is.null(n.events)) {
        ## superiority case
        if (type %in% c("superiority", "sup")) {
            power <- 1 - pnorm(logHR * sqrt(n.events * p.factor) - z_alpha)
            return(power)
        }
        ## non-inferiority case
        if (type %in% c("non-inferiority", "noninf")) {
            stopifnot(!is.null(non.inf.margin))
            power <- 1 - pnorm((logHR - non.inf.margin) * sqrt(n.events * p.factor) - z_alpha)
            return(power)
        }
    }
    stop("Wrong distribution")
}
