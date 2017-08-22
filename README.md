PhylogeneticEM
===============
[![Travis-CI Build Status](https://travis-ci.org/pbastide/PhylogeneticEM.svg?branch=develop)](https://travis-ci.org/pbastide/PhylogeneticEM)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](http://pbastide.github.io/PhylogeneticEM/)
[![codecov](https://codecov.io/gh/pbastide/PhylogeneticEM/branch/develop/graph/badge.svg)](https://codecov.io/gh/pbastide/PhylogeneticEM)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PhylogeneticEM)](https://cran.rstudio.com/web/packages/PhylogeneticEM/)

Implementation of the EM algorithm for the detection of shifts in a phylogeny.

## Installation
Stable version on the CRAN.

To get the latest (and possibly unstable) version, you can use the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools")
devtools::install_github(repo = "pbastide/PhylogeneticEM", build_vignettes = TRUE)
```

## Documentation

See package documentation (references and vignettes) here: http://pbastide.github.io/PhylogeneticEM/

(Built with [`pkgdown`](https://github.com/hadley/pkgdown)).

## References:

Bastide, P., Mariadassou, M. and Robin, S. (2017), Detection of adaptive shifts on phylogenies by using shifted stochastic processes on a tree. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4):1067-1093, [doi:10.1111/rssb.12206](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12206/abstract)

Bastide, P., Ané, C., Robin, S. and Mariadassou, M. (2017), Inference of Adaptive Shifts for Multivariate Correlated Traits. *bioRxiv preprint*, [doi:10.1101/146191](http://www.biorxiv.org/content/early/2017/06/05/146191)
