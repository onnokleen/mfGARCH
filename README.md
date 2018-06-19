[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mfGARCH)](https://cran.r-project.org/package=mfGARCH) 
[![Travis-CI Build Status](https://travis-ci.org/onnokleen/mfGARCH.svg?branch=master)](https://travis-ci.org/onnokleen/mfGARCH)
[![Coverage Status](https://img.shields.io/coveralls/onnokleen/mfGARCH.svg)](https://coveralls.io/r/onnokleen/mfGARCH?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](http://cranlogs.r-pkg.org/badges/mfGARCH)]()
# mfGARCH - mixed-frequency GARCH models

An R package for estimating multiplicative mixed-frequency GARCH models (GARCH-MIDAS) as proposed in Engle et al. (2013).

## Highlights
- A comprehensive toolbox for estimating and forecasting using GARCH-MIDAS models
- Easy to use
- Built for handling irregularly spaced mixed-frequency data

## Installation
CRAN:
```r
install.packages("mfGARCH")
```
Development version:
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
