[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mfGARCH)](https://cran.r-project.org/package=mfGARCH)
[![Travis-CI Build Status](https://travis-ci.org/onnokleen/mfGARCH.svg?branch=master)](https://travis-ci.org/onnokleen/mfGARCH)
[![Coverage Status](https://img.shields.io/coveralls/onnokleen/mfGARCH.svg)](https://coveralls.io/r/onnokleen/mfGARCH?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# mfGARCH - mixed-frequency GARCH models

An R package for estimating multiplicative mixed-frequency GARCH models (GARCH-MIDAS) as proposed in Engle et al. (2013) which is currently under development.

## Highlights
- A comprehensive toolbox for estimating, visualizing and forecasting when using GARCH-MIDAS models
- Easy to use
- Built for handling irregularly spaced mixed-frequency data

## Installation
Right now, you need the Rcpp-package to compile it on your own computer.
Development version (GitHub):
```r
# install.packages("devtools")
library(devtools)
install_github("onnokleen/mfGARCH")
```
## ToDo and Roadmap

ToDo
* Introduce non-full sample estimation by introducting sample-begin/end variable for forecasting
* Generate table function for list of mfGARCH models
* Add examples in vignette and readme

Roadmap
* 2017-06 Submission to CRAN
* 2017-06 Feature-completeness, ToDo done

Future
* Introduce other weighting schemes
* Employ more than one variable
* Broom-package plug-in
* Allow for different releases for improving real-time forecasting

## History
- January 2017: Initial submit
