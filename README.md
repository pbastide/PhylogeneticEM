PhylogeneticEM
===============
[![R-CMD-check](https://github.com/pbastide/PhylogeneticEM/workflows/R-CMD-check/badge.svg)](https://github.com/pbastide/PhylogeneticEM/actions)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](http://pbastide.github.io/PhylogeneticEM/)
[![codecov](https://codecov.io/gh/pbastide/PhylogeneticEM/branch/master/graph/badge.svg)](https://codecov.io/gh/pbastide/PhylogeneticEM)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PhylogeneticEM)](https://CRAN.R-project.org/package=PhylogeneticEM)

Implementation of the EM algorithm for the detection of shifts in a phylogeny.

## Installation
Stable version on the [CRAN](https://cran.rstudio.com/web/packages/PhylogeneticEM/).

To get the latest (and possibly unstable) version, you can use the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools")
devtools::install_github(repo = "pbastide/PhylogeneticEM", ref = "develop")
```

## Documentation

See package documentation (references and vignettes) here: http://pbastide.github.io/PhylogeneticEM/

(Built with [`pkgdown`](https://github.com/hadley/pkgdown)).

## References:

Bastide, P., Mariadassou, M. and Robin, S. (2017), Detection of adaptive shifts on phylogenies by using shifted stochastic processes on a tree. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4):1067-1093, [doi:10.1111/rssb.12206](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12206/abstract).

Bastide, P., Ané, C., Robin, S. and Mariadassou, M. (2018), Inference of Adaptive Shifts for Multivariate Correlated Traits. *Systematic Biology*, 67(4), 662–680. [doi:10.1093/sysbio/syy005](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy005/4827615?guestAccessKey=fba26a20-0579-4721-ad76-8e669489539a).
