# v2.2.1
  - added option `type = "paired"` to `power.z.test`
  - Changes in `confIntIndependentAUCDiff`:
    - fixed error where setting `type = "Wald"` would return the same result as `type = "logit"` due to an issue with uppercase and lowercase characters.
    - changed choices for the `type` argument to lowercase `"wald"` and `"logit"`. However, the uppercase variants still work.
    - removed code repetitions in function.
  - change dependency structure such that the package only depends on R itself. 
    - All other dependencies are moved from `Depends` to the `Imports` section of the `DESCRIPTION` file. This avoids potential issues when other packages mask functions that `biostatUZH` relies on.
    - Calls to functions from other packages are now clearly labelled through consequent use of the `::` operator.

# v2.2
  - added function `power.z.test` which is based on `power.t.test`.

# v2.1
  - The package now reexports the functions `ci2p`, `ci2estimate`, `ci2se`, `ci2z`, `z2p`, `p2z` from ReplicationSuccess

# v2.0.4
  - transfer ownership from https://github.com/felix-hof/biostatUZH to https://github.com/EBPI-Biostatistics/biostatUZH
  - move `prodlim` package from `Suggests` to `Imports` section

# v2.0.3
  - allow NAs in formatPval() function

# v2.0.2
  - update links
  - fix errors in tableRegression() input checks

# v2.0.1
  - migrate package from https://github.com/florafauna/biostatUZH to https://github.com/felix-hof/biostatUZH

# v2.0.0
  - add NEWS.md file
  - add Makefile to support building and testing the package
  - migrate from R-forge to github repo:
    https://github.com/florafauna/biostatUZH
  - improve R coding to make it more stable and readable
  - unify function and argument names:
    function names start with lower case and use camelStyle,
    argument names are lower case
  - improve documentation and migrate it to roxygen2 format
  - implement input check for all functions
  - combine related functions into one R file
  - remove unused function as discussed with Leo


# v1.8.0
  - Available from R-forge:
    https://r-forge.r-project.org/R/?group_id=2036

