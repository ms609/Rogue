# Rogue

[![codecov](https://codecov.io/gh/ms609/TreeSearch/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/RogueNaRok)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/RogueNaRok)](https://cran.r-project.org/package=RogueNaRok)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/RogueNaRok)](https://cran.r-project.org/package=RogueNaRok)
[![DOI](https://zenodo.org/badge/98171XXX642.svg)](https://zenodo.org/badge/latestdoi/98171642)[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)

This package provides an R interface to [RogueNaRok](https://rnr.h-its.org/about),
an algorithm to prune rogue taxa and improve the resolution of phylogenetic
consensus trees (Aberer _et al._ 2013, [doi:10.1093/sysbio/sys078](https://dx.doi.org/10.1093/sysbio/sys078)).


# Installation

<!--Install and load the stable version from CRAN, and launch the GUI, as follows:-->
Install and load the stable version from CRAN as follows:
```r
install.packages('Rogue')
library('Rogue')
```

Install the develoment version from GitHub as follows:
```r
devtools::install_github('ms609/Rogue')
library('Rogue')
```

# References

A.J. Aberer, D. Krompass, A. Stamatakis (2017): Pruning Rogue Taxa Improves
  Phylogenetic Accuracy: An Efficient Algorithm and Webservice, Systematic Biology 62(1):
  162-166, doi:[10.1093/sysbio/sys078](https://dx.doi.org/10.1093/sysbio/sys078).
