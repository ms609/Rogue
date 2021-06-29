# Rogue

 [![Codecov test coverage](https://codecov.io/gh/ms609/Roguer/branch/main/graph/badge.svg)](https://codecov.io/gh/ms609/Roguer?branch=main)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/Rogue)](https://cran.r-project.org/package=Rogue)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/Rogue)](https://cran.r-project.org/package=Rogue)
[![DOI](https://zenodo.org/badge/376830950.svg)](https://zenodo.org/badge/latestdoi/376830950)[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)

An R interface to [RogueNaRok](https://rnr.h-its.org/about),
an algorithm to prune rogue taxa and improve the resolution of phylogenetic
consensus trees (Aberer _et al._ 2017).

# Installation

<!--Install and load the stable version from CRAN, and launch the GUI, as follows:-->
Install and load the stable version from CRAN as follows:
```r
install.packages('Rogue')
library('Rogue')
```

Install the development version from GitHub as follows:
```r
devtools::install_github('ms609/Rogue')
library('Rogue')
```

If you find the package useful in your work, please consider citing 
Aberer _et al._ (2017) (description and implementation of original algorithm)
and Smith (2021) (R interface to Aberer's C code).

# References

A.J. Aberer, D. Krompass, A. Stamatakis (2017): Pruning Rogue Taxa Improves
  Phylogenetic Accuracy: An Efficient Algorithm and Webservice, Systematic Biology 62(1):
  162-166, doi:[10.1093/sysbio/sys078](https://dx.doi.org/10.1093/sysbio/sys078).

M.R. Smith (2021): Rogue: Identify Rogue Taxa in Sets of Phylogenetic Trees.
  Zenodo,
  doi:[10.5281/zenodo.5037327](https://dx.doi.org/10.5281/zenodo.5037327).
