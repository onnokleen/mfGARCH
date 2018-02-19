[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mfGARCH)](https://cran.r-project.org/package=mfGARCH) 
[![Travis-CI Build Status](https://travis-ci.org/onnokleen/mfGARCH.svg?branch=master)](https://travis-ci.org/onnokleen/mfGARCH)
[![Coverage Status](https://img.shields.io/coveralls/onnokleen/mfGARCH.svg)](https://coveralls.io/r/onnokleen/mfGARCH?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# mfGARCH - mixed-frequency GARCH models

An R package for estimating multiplicative mixed-frequency GARCH models (GARCH-MIDAS) as proposed in Engle et al. (2013).

## Highlights
- A comprehensive toolbox for estimating and forecasting using GARCH-MIDAS models
- Easy to use
- Built for handling irregularly spaced mixed-frequency data

## Installation
Right now, you need the Rcpp-package to compile it on your own computer.
Development version (GitHub):
```r
# Install package via devtools
# install.packages("devtools")
library(devtools)
install_github("onnokleen/mfGARCH")

# Example
library(mfGARCH)
# df_financial
fit_mfgarch(data = df_financial, y = "return", x = "nfci", low.freq = "week", K = 52)
```
## ToDo and Roadmap

ToDo
* Add examples in vignette and readme
* Add generic plot function

Roadmap
* 2018-02 Submission to CRAN

Future
* Broom-package plug-in

## History
- January 2017: Initial submit
