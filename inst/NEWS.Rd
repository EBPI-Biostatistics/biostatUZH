%% Pre-defined commands: \R, \code, \acronym, \url, \file, \pkg
%% Since R 3.2.0, there are some additional system Rd macros available,
%% e.g., \CRANpkg and \doi. See the definitions in the file
%% file.path(R.home("share"), "Rd", "macros", "system.Rd")

\name{NEWS}
\title{News for Package 'biostatUZH'}

\section{Changes in biostatUZH version 2.2.2 (2022-11-04)}{
  \itemize{
    \item Made \code{behrens.test} a generic with a \code{default} method and a method for \code{formula} objects.
    \item Added function \code{formatPvalStrict} which allows the user to fix the number of digits after the decimal point.
    \item Added an argument \code{strict} to the function \code{tableRegression} which indicates whether *p*-values should be formatted with \code{formatPvalStrict} or \code{formatPval}.
  }
}

\section{Changes in biostatUZH version 2.2.1 (2022-11-01)}{
  \itemize{
    \item Renamed function \code{behrensTest} to \code{behrens.test} and return an object of class \code{htest} such that the print methods from the \pkg{stats}-package can be used on it.
  }
}

\section{Changes in biostatUZH version 2.2 (2022-10-24)}{
  \itemize{
    \item Added the function \code{power.z.test} that allows sample size calculations based on normal approximations.
  }
}

\section{Changes in biostatUZH version 2.1 (2022-09-28)}{
  \itemize{
    \item The package now reexports the functions \code{ci2p}, \code{ci2estimate}, \code{ci2se}, \code{ci2z}, \code{z2p}, \code{p2z} from ReplicationSuccess
  }
}

\section{Changes in biostatUZH version 2.0.4 (2022-09-14)}{
  \itemize{
    \item transfer package from \url{https://github.com/felix-hof/biostatUZH} to \url{https://github.com/EBPI-Biostatistics/biostatUZH}
    \item move \pkg{prodlim} package from \code{Suggests} to \code{Imports} section in the \file{DESCRIPTION} file. This fixes an error about the missing package in the installation process.
  }
}

\section{Changes in biostatUZH version 2.0.3 (2022-03-07)}{
  \itemize{
    \item allow NAs in formatPval() function
  }
}

\section{Changes in biostatUZH version 2.0.2 (2021-10-27)}{
  \itemize{
    \item update links in documentation
    \item fix errors in tableRegression() input checks
    \item restore commit history from R-forge
  }
}

\section{Changes in biostatUZH version 2.0.1 (2021-10-07)}{
  \itemize{
    \item migrate package from \url{https://github.com/florafauna/biostatUZH} to \url{https://github.com/felix-hof/biostatUZH}
  }
}

\section{Changes in biostatUZH version 2.0.0 (2021-09-14)}{
  \itemize{
    \item The function \code{tableOR()} has been added to create
    results tables from logistic regression fits for LaTeX.
    \item add NEWS.md file
    \item add Makefile to support bilding and testing the package
    \item migrate from R-forge to github repo: \url{https://github.com/florafauna/biostatUZH}
    \item improve R coding to make it more stable and readable
    \item unify function and argument names:
      \itemize{
        \item function names start with lower case and use camelStyle,
        \item argument names are lower case
      }
    \item improve documentation and migrate it to roxygen2 format
    \item implement input check for all functions
    \item combine related functions into one R file
    \item remove unused function as discussed with Leo
  }
}

\section{Changes in biostatUZH version 1.8.0 (2018-10-18)}{
  \itemize{
    \item The function \code{tableOR()} has been added to create
    results tables from logistic regression fits for LaTeX.
  }
}

\section{Changes in biostatUZH version 1.7.0 (2016-11-18)}{
  \itemize{
    \item The function \code{confIntAUC()} has been replaced by a new
    version (with a different result) using the new function
    \code{standardErrorAUC()}. The old implementation of
    \code{confIntAUC()} is available as \code{summaryROC()}.
    \item Several CI's for operating characteristics of diagnostic tests
    have been implemented: \code{confIntDiagnostic()},
    \code{confIntIndependentDiagnostic()},
    \code{confIntPairedDiagnostic()}.
    Furthermore, \code{confIntRiskRatio()} has been added.
    \item The new function \code{behrensTest()} computes Behrens'
    version of the t-test for unequal variances.
    \item \code{formatPval()} now supports missing values.
  }
}

\section{Changes in biostatUZH version 1.6.0 (2016-02-24)}{
  \itemize{
    \item New function \code{jeffreys()} supersedes \code{jeffrey()}
    and returns the posterior median as the point estimate instead of
    \eqn{x/n}. Furthermore, the CI limits for the special cases
    \eqn{x=0} and \eqn{x=n} are no longer forced to 0 and 1,
    respectively.
    \item \code{wald()} now truncates confidence limits at 0 and 1.
    \item Functions \code{confIntProportion()}, \code{wald()},
    \code{wilson()}, \code{agresti()}, and \code{clopperPearson()} now
    check their input arguments.
    \item In \code{confIntProportion()$CIs$type}, the name
    \code{"clopperPearson"} has been changed to \code{"ClopperPearson"}.
    \item New alias \code{expit()} for \code{ilogit()}.
  }
}

\section{Changes in biostatUZH version 1.5.3 (2015-10-16)}{
  \itemize{
    \item Added a convenient shortcut \code{plotcorr(..., col = TRUE)}
    to fill the ellipses with colors mapped to the correlations (which
    is the new default).
  }
}

\section{Changes in biostatUZH version 1.5.2 (2015-09-09)}{
  \itemize{
    \item Fixed yet another \code{conf.level} bug
    (sitting in \code{WeibullReg()} since version 1.1.0).
  }
}

\section{Changes in biostatUZH version 1.5.1 (2015-05-06)}{
  \itemize{
    \item \code{jeffrey()} and thus \code{confIntProportion()} now also
    work in the boundary cases \code{x = 0} and \code{x = n}.
  }
}

\section{Changes in biostatUZH version 1.5.0 (2015-02-12)}{
  \itemize{
    \item New function \code{formatCI()}, a re-implementation of
    \code{displayCI()} from package \pkg{reporttools}.
    \item New function \code{plotcorr()}, a fork from package
    \pkg{ellipse}, which enables ellipses in the lower and numbers in
    the upper triangle (or vice versa).
    \item Bugfix in \code{confIntAUC()}: the argument \code{conf.level}
    was ignored when calculating Wilson confidence intervals for
    specificity, sensitivity, PPV, and NPV.
    \item \code{tableRegression()} ...
    \itemize{
      \item supports negative binomial models fit by \code{glm.nb} from
      package \pkg{MASS}.
      \item passes additional arguments to \code{print.xtable} and
      sanitizes row names by default.
    }
    \item Migrated package development to the R-Forge project \dQuote{ebuzh}.
    \item Accommodated R CMD check issues.
    \item Updated institute name.
    \item Standardized \file{NEWS} file format.
    \item The \pkg{boot} and \pkg{reporttools} packages are no longer
    automatically attached when attaching package \pkg{biostatUZH}.
  }
}

\section{Changes in biostatUZH version 1.4.9 (2014-04-15)}{
  \itemize{
    \item Bugfix in \code{confIntIndependentProportion()}: the Newcombe interval
    was wrong due to wrong internal use of the \code{conf.level} argument
    (spotted by Andrea Kraus). The bug was probably introduced in
    version 1.1.0, when all \code{alpha} arguments were changed to
    \code{conf.level}.
  }
}

\section{Changes in biostatUZH version 1.4.8 (2013-11-21)}{
  \itemize{
    \item improved code and doc of \code{formatPval()} and \code{logit()}
    \item fixed \file{NAMESPACE} issues
    \item moved most packages from Depends and Imports to Suggests
    \item consistently use upper-case \file{.R} as the source file extension
  }
}

\section{Changes in biostatUZH version 1.4.7 (2013-04-08)}{
  \itemize{
    \item \code{tableRegression()}: bugfix (related to \code{intercept} argument)
  }
}

\section{Changes in biostatUZH version 1.4.6 (2013-04-05)}{
  \itemize{
    \item \code{tableRegression()}: removed default
    \code{"coefficients"} from all non-lm objects
  }
}

\section{Changes in biostatUZH version 1.4.5 (2013-04-04)}{
  \itemize{
    \item \code{tableRegression()}: completely updated
  }
}

\section{Changes in biostatUZH version 1.4.4 (2013-03-07)}{
  \itemize{
    \item \code{tableRegression()}: option \code{"text"} added (before
    there was only \code{"text.ci"})
  }
}

\section{Changes in biostatUZH version 1.4.3 (2013-01-18)}{
  \itemize{
    \item minor bug fixes
  }
}

\section{Changes in biostatUZH version 1.4.2 (2012-11-26)}{
  \itemize{
    \item \code{tableRegression()}: removed automatic labeling of HR, OR and RR
  }
}

\section{Changes in biostatUZH version 1.4.1 (2012-11-07)}{
  \itemize{
    \item \code{tableRegression()}: added argument \code{biostatUZH.style}
  }
}

\section{Changes in biostatUZH version 1.4.0 (2012-10-09)}{
  \itemize{
    \item new functions \code{tableRegression()} and \code{formatPval()}
  }
}

\section{Changes in biostatUZH version 1.3.0 (2012-10-08)}{
  \itemize{
    \item new function \code{confIntRiskDiff()}
  }
}

\section{Changes in biostatUZH version 1.2.0 (2012-09-21)}{
  \itemize{
    \item \code{confIntAUCbinorm()} has been extended with Wald CIs
  }
}

\section{Changes in biostatUZH version 1.1.2 (2011-04-28)}{
  \itemize{
    \item corrected minor bug in \code{quantileKM()} (thanks to SH)
  }
}

\section{Changes in biostatUZH version 1.1.1 (2011-01-13)}{
  \itemize{
    \item Added function \code{sampleSizePhase2BinaryTest()}
  }
}

\section{Changes in biostatUZH version 1.1.0 (2010-11-23)}{
  \itemize{
    \item Added functions \code{formatPercent()} and
    \code{confIntIndependenProportion()} (thanks to LH)
    \item changed all arguments \code{alpha=0.05} to
    \code{conf.level=0.95} in functions for confidence intervals
  }
}

\section{Changes in biostatUZH version 1.0.0 (2010-06-25)}{
  \itemize{
    \item renamed package to \pkg{biostatUZH}
    \item added functions \code{confIntProportion()} (collects all
    confidence intervals for a single proportion),
    \code{sampleSizeMcNemar()}, \code{confIntPairedProportion()},
    \code{confIntMedian()}, and \code{quantileCumInc()}
    \item fixed bug in \code{confIntKappa()} for \code{"bca"} confidence
    interval (remark SH)
    \item \code{?wilson}, \code{?wald}, etc now work (remark by LH)
  }
}

\section{Changes in biostatUZH version 0.0.0.9000 (2010-03-17)}{
  \itemize{
    \item 17.03.2010: added function \code{survrep()} of SH
    \item 03.06.2009: integrated Weibull-functions of SH
    \item 15.05.2009: added \code{confIntFisherTrafo()}
  }
}
