# v2.2.6
  - Renamed argument `alpha` of function `survEvent()` to `sig.level`. This fixes a small inconsistency in argument naming.

# v2.2.5
  - Added function confIntRateRatio.

# v2.2.4
  - Minor fix: return both, *p*-value and test statistic in `gailSimon()`.

# v2.2.3
  - Added implementation of the Gail-Simon test for qualitative interaction in the function `gailSimon()`.

# v2.2.2
  - Made `behrens.test()` a generic with a `default` method and a method for `formula` objects.
  - Added function `formatPvalStrict()` which allows the user to fix the number of digits after
  the decimal point.
  - Added an argument `strict` to the function `tableRegression()` which indicates whether *p*-values should be formatted with `formatPvalStrict()` or `formatPval()`.

# v2.2.1
  - Renamed function `behrensTest()` to `behrens.test()` and return an object of class `htest` such that the print methods from the `stats`-package can be used on it.

# v2.2
  - Added function `power.z.test()` which is based on `power.t.test()`.

# v2.1
  - The package now reexports the functions `ci2p()`, `ci2estimate()`, `ci2se()`, `ci2z()`, `z2p()`, `p2z()` from ReplicationSuccess

# v2.0.4
  - Transfer ownership from https://github.com/felix-hof/biostatUZH to https://github.com/EBPI-Biostatistics/biostatUZH
  - Move `prodlim` package from `Suggests` to `Imports` section

# v2.0.3
  - Allow `NA`s in `formatPval()` function

# v2.0.2
  - Update links
  - Fix errors in `tableRegression()` input checks

# v2.0.1
  - Migrate package from https://github.com/florafauna/biostatUZH to https://github.com/felix-hof/biostatUZH

# v2.0.0
  - Add NEWS.md file
  - Add Makefile to support bilding and testing the package
  - Migrate from R-forge to github repo:
    https://github.com/florafauna/biostatUZH
  - Improve R coding to make it more stable and readable
  - Unify function and argument names:
    function names start with lower case and use camelStyle,
    argument names are lower case
  - Improve documentation and migrate it to roxygen2 format
  - Implement input check for all functions
  - Combine related functions into one R file
  - Remove unused function as discussed with Leo


# v1.8.0
  - Available from R-forge:
    https://r-forge.r-project.org/R/?group_id=2036

