# Rogue

 [![Codecov test coverage](https://codecov.io/gh/ms609/Rogue/branch/main/graph/badge.svg)](https://codecov.io/gh/ms609/Rogue?branch=main)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/Rogue)](https://cran.r-project.org/package=Rogue)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/Rogue)](https://cran.r-project.org/package=Rogue)
[![DOI](https://zenodo.org/badge/376830950.svg)](https://zenodo.org/badge/latestdoi/376830950)
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)

"Rogue" implements approaches to identify rogue taxa in phylogenetic analysis.
Rogues are wildcard leaves whose uncertain position reduces the resolution of
consensus trees. Consensus trees that omit rogue taxa can be more informative.

"Rogue" allows the user to select a concept of "information" by which the
quality of consensus trees should be evaluated, and a heuristic approach
by which rogue taxa should be identified.

Rogue detection using the phylogenetic and clustering information content
measures (Smith, forthcoming) is implemented using a quick heuristic that drops
the least "stable" leaves one at a time,
using an _ad hoc_ definition of stability (Smith, forthcoming);
and by a more exhaustive (and time-consuming) approach that considers dropping
all possible sets of up to _n_ leaves (Aberer _et al._ 2013).

The latter heuristic is implemented for the relative bipartition 
"information" content and Pattengale's criterion
_via_ [RogueNaRok](https://rnr.h-its.org/about) (Aberer _et al._ 2013).


# Installation

Install and load the stable version from CRAN as normal:
```r
install.packages('Rogue')
library('Rogue')
```

Implementations of the phylogenetic and clustering information criteria are
not yet available on CRAN.

Install the development version from GitHub with 
`devtools::install_github('ms609/Rogue')`.


# Citing 'Rogue'

If you find this package useful in your work, Please consider citing
Smith (2021).

To cite the underlying methods, please cite Aberer _et al._ (2013) ('RogueNaRok')
or Smith (forthcoming), as appropriate.


# References

A.J. Aberer, D. Krompass, A. Stamatakis (2013): Pruning Rogue Taxa Improves
  Phylogenetic Accuracy: An Efficient Algorithm and Webservice, Systematic Biology 62(1):
  162-166, doi:[10.1093/sysbio/sys078](https://dx.doi.org/10.1093/sysbio/sys078).

M.R. Smith (2021): Rogue: Identify Rogue Taxa in Sets of Phylogenetic Trees.
  Zenodo,
  doi:[10.5281/zenodo.5037327](https://dx.doi.org/10.5281/zenodo.5037327).

M.R. Smith (forthcoming): Improving consensus trees by detecting rogue taxa.
