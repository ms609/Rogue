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
measures (Smith, 2022) is implemented using a quick heuristic that drops
the least "stable" leaves one at a time,
using an _ad hoc_ definition of stability (Smith, 2022);
and by a more exhaustive (and time-consuming) approach that considers dropping
all possible sets of up to _n_ leaves (Aberer _et al._ 2013).

The latter heuristic is implemented for the relative bipartition 
"information" content and Pattengale's criterion
_via_ [RogueNaRok](https://rnr.h-its.org/about) (Aberer _et al._ 2013).


[![Detecting rogue taxa with information theory](man/figures/Rogue_talk.png)](https://durham.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=86c175f1-6e20-499c-bcf2-adeb0137a4a7)

# Installation

Install and load the stable version from CRAN as normal:
```r
install.packages("Rogue")
library("Rogue")
```

Install the development version from GitHub with 
`devtools::install_github("ms609/Rogue", args="--recursive")`.
(Requires [git](https://git-scm.com/) to be installed and added to
your PATH system environment variable; you may also require the "curl" R package.)


# Citing 'Rogue'

If you find this package useful in your work, please consider citing
Smith (2021).

To cite the underlying methods, please cite Aberer _et al._ (2013) ('RogueNaRok')
or Smith (2022), as appropriate.


# References

A.J. Aberer, D. Krompass, A. Stamatakis (2013): Pruning rogue taxa improves
  phylogenetic accuracy: an efficient algorithm and webservice. _Systematic Biology_ 62(1):
  162-166, [doi:10.1093/sysbio/sys078](https://dx.doi.org/10.1093/sysbio/sys078).

M.R. Smith (2021): Rogue: Identify rogue taxa in sets of phylogenetic trees.
  _Zenodo_,
  [doi:10.5281/zenodo.5037327](https://dx.doi.org/10.5281/zenodo.5037327).

M.R. Smith (2022): Using information theory to detect rogue taxa and improve
  consensus trees. _Systematic Biology_, syab099,
  [doi:10.1093/sysbio/syab099](https://dx.doi.org/10.1093/sysbio/syab099)
