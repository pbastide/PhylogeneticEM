## Test environments
* local: macOS 10.15, R 4.2.1
* GitHub Actions:
  * macOS-latest: release
  * windows-latest: release
  * ubuntu-latest: devel, release and oldrel
* win-builder: devel and release

## R CMD check results
There were no ERRORs or WARNINGs.
There were two NOTEs:
```
checking installed package size ... NOTE
  installed size is 5.0Mb
  sub-directories of 1Mb or more:
    libs 3.5Mb
```
This package uses `Rcpp` and `RcppArmadillo`.

```
* checking CRAN incoming feasibility ... [18s] NOTE
Maintainer: 'Paul Bastide <paul.bastide@m4x.org>'

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1111/rssb.12206
    From: inst/CITATION
    Status: 503
    Message: Service Unavailable
```
This doi and this address should be valid.

## Downstream dependencies
There are currently no downstream dependencies for this package.