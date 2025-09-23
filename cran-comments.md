## Test environments
* local: macOS 14.7.8, R 4.5.1
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
  installed size is 5.1Mb
  sub-directories of 1Mb or more:
    libs 3.5Mb
```
This package uses `Rcpp` and `RcppArmadillo`.

## Downstream dependencies
There are currently no downstream dependencies for this package.
