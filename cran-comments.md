## Test environments
* local: macOS 10.14, R 4.0.3
* GitHub Actions:
  * macOS-latest: release
  * windows-latest: release
  * ubuntu-latest: devel, release and oldrel
* win-builder: devel and release

## R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE:
```
checking installed package size ... NOTE
  installed size is 5.0Mb
  sub-directories of 1Mb or more:
    libs 3.5Mb
```
This package uses `Rcpp` and `RcppArmadillo`.

## Downstream dependencies
There are currently no downstream dependencies for this package.