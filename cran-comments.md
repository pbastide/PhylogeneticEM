## Test environments
* macOS 10.12 (local), R 3.5.1
* ubuntu 16.04 (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs.

There was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Paul Bastide <paul.bastide@m4x.org>'

Possibly mis-spelled words in DESCRIPTION:
  Bastide (11:20, 12:5)
  al (11:31, 12:16)
  et (11:28, 12:13)
  
Those words are correctly spelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Comments on resubmission
The checking and vignette-building time was reduced by half, to about 250s (40 + 40 + 170, on win-builder with R-devel).
Two references were added to the Description field.
Thank you !