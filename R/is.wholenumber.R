
# Function taken from the R base manual (src/library/base/man/integer.Rd)
# Copyright 1995-2010 R Core Team
# Distributed under GPL 2 or later

is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
    abs(x - round(x)) < tol
