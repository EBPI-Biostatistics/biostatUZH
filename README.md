# biostatUZH

Provides R functions developed at the Department of Biostatistics,
Epidemiology, Biostatistics and Prevention Institute,
University of Zurich, Zurich, Switzerland.


Currently implemented topics:
* several CI's for proportions,
* several CI's for diagnostic tests,
* Fagan-Nomogram,
* bootstrap CI's for the kappa coefficient,
* intraclass-correlation coefficients including CI's,
* agreement for continuous measurements (Bland-Altman plot),
* CI for the area under the curve (AUC),
* Mantel-Haenszel estimator,
* Mantel-Cox hazard ratio estimator,
* CI for the Kaplan-Meier estimate at given time points,
* CI for quantile of a Kaplan-Meier or cumulative incidence estimate,
* sample size computations for two-sample Mann-Whitney test,
* sample size computations for the McNemar test,
* sample size computations for survival outcomes, 
* binary diagnostic test (via normal approximation and simulation),
* natural re-parametrization of Weibull output from survreg,
* hazard ratio and event time ratio interpretations,
* a plot to check the adequacy of the Weibull model.

## Installation

```r
devtools::install_github(repo = "florafauna/biostatUZH")
```

## Vignettes

```r
vignette("weibull"
vignette("aucbinormal")
```

## Website

* http://www.biostat.uzh.ch
