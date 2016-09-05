\name{biostatUZH-package}
\alias{biostatUZH-package}
\alias{biostatUZH}
\docType{package}
\keyword{package}

\title{Misc Tools of the Department of Biostatistics, EBPI, University of Zurich}

\description{
  The package \pkg{biostatUZH} collects various functions developed at
  the Department of Biostatistics at the Epidemiology, Biostatisics and
  Prevention Institute, University of Zurich, Switzerland.
  Currently implemented topics:
  several CI's for proportions,
  several CI's for operating characteristics of diagnostic tests,
  Fagan-Nomogram,
  bootstrap CI's for the kappa coefficient,
  intraclass-correlation coefficients including CI's,
  agreement for continuous measurements (Bland-Altman plot),
  CI for the area under the curve (AUC),
  Mantel-Haenszel estimator,
  Mantel-Cox hazard ratio estimator,
  CI for the Kaplan-Meier estimate at given time points,
  CI for quantile of a Kaplan-Meier or cumulative incidence estimate,
  sample size computations for two-sample Mann-Whitney test,
  McNemar test,
  a binary diagnostic test (via normal approximation and simulation),
  natural re-parametrization of Weibull output from survreg,
  hazard ratio and event time ratio interpretations,
  plot to check the adequacy of the Weibull model.
}

\details{
\tabular{ll}{
Package:  \tab biostatUZH \cr
URL:      \tab \url{http://www.biostat.uzh.ch} \cr
License:  \tab GPL (>=2)  \cr
}

Functions that implement methods that are not implemented in \R
elsewhere are, to the best of our knowledge: \cr
\code{\link{BlandAltman}} \cr
\code{\link{computeICCrater}} \cr
\code{\link{confIntAUC}} \cr
\code{\link{confIntAUCbinorm}} \cr
\code{\link{confIntDiagnostic}} \cr
\code{\link{confIntICC}} \cr
\code{\link{confIntKappa}} \cr
\code{\link{confIntPairedDiagnostic}} \cr
\code{\link{confIntPairedProportion}} \cr
\code{\link{confIntIndependentDiagnostic}} \cr
\code{\link{confIntIndependentProportion}} \cr
\code{\link{confIntRiskDiff}} \cr
\code{\link{ConvertWeibull}} \cr
\code{\link{faganPlot}} \cr
\code{\link{populationSamplePlot}} \cr
\code{\link{quantileCumInc}} \cr
\code{\link{sampleSizeMcNemar}} \cr
\code{\link{WeibullReg}} \cr
\code{\link{WeibullDiag}} \cr
\code{\link{sampleSizePhase2BinaryTest}} \cr
\code{\link{sampleSizePhase2BinaryTestSimulated}}
}

\author{
  Main authors: Leonhard Held, Kaspar Rufibach
  
  Secondary authors: Sarah Haile, Sebastian Meyer, Sina Rueeger, 

  Further contributions by: Andrea Riebler, Daniel Sabanes Bove

  Maintainer: Leonhard Held \email{leonhard.held@uzh.ch}
}
