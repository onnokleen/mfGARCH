[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mfGARCH)](https://cran.r-project.org/package=mfGARCH) 
[![Travis-CI Build Status](https://travis-ci.org/onnokleen/mfGARCH.svg?branch=master)](https://travis-ci.org/onnokleen/mfGARCH)
[![Coverage Status](https://img.shields.io/coveralls/onnokleen/mfGARCH.svg)](https://coveralls.io/r/onnokleen/mfGARCH?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://cranlogs.r-pkg.org/badges/mfGARCH)](https://cranlogs.r-pkg.org/badges/mfGARCH)
# mfGARCH - mixed-frequency GARCH models

An R package for estimating GARCH-MIDAS (MIxed-DAta-Sampling) models (Engle, Ghysels and Sohn, 2013, [doi:10.1162/REST_a_00300](https://doi.org/10.1162/REST_a_00300)) and related statistical inference, accompanying the paper "Two are better than one: Volatility forecasting using multiplicative component GARCH models" by Conrad and Kleen (2020, [doi:10.1002/jae.2742](https://doi.org/10.1002/jae.2742)). The GARCH-MIDAS model decomposes the conditional variance of (daily) stock returns into a short- and long-term component, where the latter may depend on an exogenous covariate sampled at a lower frequency.

## Highlights
- A comprehensive toolbox for estimating and forecasting using GARCH-MIDAS models
- Easy to use, both with one or two explanatory covariates
- Built for handling irregularly spaced mixed-frequency data

Please cite as

> Conrad, Christian and Kleen, Onno (2020). Two are better than one: Volatility forecasting using multiplicative component GARCH-MIDAS models. Journal of Applied Econometrics 35: 19&ndash;45.

and

> Kleen, Onno (2020). mfGARCH: Mixed-Frequency GARCH Models. R package version 0.2.0.

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
```

## Example
```r
library(mfGARCH)
# df_financial
fit_mfgarch(data = df_financial, y = "return", x = "nfci", low.freq = "week", K = 52)
```
