%% Pre-defined commands: \R, \code, \acronym, \url, \file, \pkg
%% Since R 3.2.0, there are some additional system Rd macros available,
%% e.g., \CRANpkg and \doi. See the definitions in the file
%% file.path(R.home("share"), "Rd", "macros", "system.Rd")

\name{NEWS}
\title{News for Package 'biostatUZH'}

\section{Changes in biostatUZH version 1.5.0 (2015-02-10)}{
  \itemize{
    \item Bugfix in \code{confIntAUC()}: the argument \code{conf.level}
    was ignored when calculating Wilson confidence intervals for
    specificity, sensitivity, PPV, and NPV.
    \item Migrated package development to the R-Forge project \dQuote{ebuzh}.
    \item Accommodated R CMD check issues.
    \item Updated institute name.
    \item Standardized \file{NEWS} file format.
    \item The \pkg{boot} package is no longer attached (but only imported).
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
