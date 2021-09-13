# biostatUZH

Provides R functions developed at the Department of Biostatistics,
Epidemiology, Biostatistics and Prevention Institute,
University of Zurich, Zurich, Switzerland.
http://www.biostat.uzh.ch


Covered topics:
* CIs for proportions,
* CIs for diagnostic tests,
* bootstrap CIs for the kappa coefficient,
* intraclass-correlation coefficients including CIs,
* agreement for continuous measurements (Bland-Altman plot),
* CI for the area under the curve (AUC),
* Mantel-Haenszel estimator,
* Mantel-Cox hazard ratio estimator,
* CI for the Kaplan-Meier estimate at given time points,
* CI for quantile of a Kaplan-Meier or cumulative incidence estimate,
* sample size computations for two-sample Mann-Whitney test,
* sample size computations for the McNemar test,
* sample size computations for survival outcomes, 
* binary diagnostic test
* natural re-parametrization of Weibull output from survreg,
* hazard ratio and event time ratio interpretations,
* a plot to check the adequacy of the Weibull model.

## Get started:

Install the package
```r
devtools::install_github(repo = "florafauna/biostatUZH")
```

Load package and list all provided functions
```r
library(biostatUZH)
ls("package:biostatUZH")
```

Get the documentation of a function
```r
?confIntProportion
```

Access vignettes
```r
vignette("weibull")
vignette("aucbinormal")
```