PhylogeneticEM
===============
[![Travis-CI Build Status](https://travis-ci.org/pbastide/PhylogeneticEM.svg?branch=master)](https://travis-ci.org/pbastide/PhylogeneticEM)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](http://pbastide.github.io/PhylogeneticEM/)
[![codecov](https://codecov.io/gh/pbastide/PhylogeneticEM/branch/master/graph/badge.svg)](https://codecov.io/gh/pbastide/PhylogeneticEM)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PhylogeneticEM)](https://cran.rstudio.com/web/packages/PhylogeneticEM/)

Implementation of the EM algorithm for the detection of shifts in a phylogeny.

## Installation
Stable version to be released on the CRAN soon.

To get the latest (and possibly unstable) version, you can use the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools")
devtools::install_github(repo = "pbastide/PhylogeneticEM", build_vignettes = TRUE)
```

## Documentation

See package documentation (references and vignettes) here: http://pbastide.github.io/PhylogeneticEM/

(Built with [`pkgdown`](https://github.com/hadley/pkgdown)).

## Reference:

Bastide, P., Mariadassou, M. and Robin, S. (2016), Detection of adaptive shifts on phylogenies by using shifted stochastic processes on a tree. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, [doi:10.1111/rssb.12206](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12206/abstract)
