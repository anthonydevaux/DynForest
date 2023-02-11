
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/DynForest)](https://CRAN.R-project.org/package=DynForest)
[![Downloads](https://cranlogs.r-pkg.org/badges/DynForest?color=blue)](https://www.r-pkg.org/pkg/DynForest)
[![R build
status](https://github.com/anthonydevaux/DynForest/workflows/R-CMD-check/badge.svg)](https://github.com/anthonydevaux/DynForest/actions)
<!-- badges: end -->

## Description

`DynForest` is a R package aiming to predict an outcome using
multivariate longitudinal predictors. The method is based on random
forest principle where the longitudinal predictors are modeled through
the random forest. `DynForest` currently supports continuous,
categorical and survival outcome. The methodology is fully described for
a survival outcome in the paper:

> Devaux A., Helmer C., Genuer R. and Proust-Lima C. (2022). Random
> survival forests with multivariate longitudinal endogenous covariates.
> *arXiv*. [\<doi:
> 10.48550/arXiv.2208.05801\>](https://doi.org/10.48550/arXiv.2208.05801)

`DynForest` user guide is also available in the paper:

> Devaux A., Proust-Lima C. and Genuer R. (2023). Random Forests for
> time-fixed and time-dependent predictors: The DynForest R package.
> *arXiv*.
> [\<doi:10.48550/arXiv.2302.02670\>](https://doi.org/10.48550/arXiv.2302.02670)

## Installation

`DynForest` package version 1.1.1 could be install from the
[CRAN](https://cran.r-project.org/package=DynForest) with:

``` r
install.packages("DynForest")
```

Development version of `DynForest` is also available from
[GitHub](https://github.com/anthonydevaux/DynForest) with:

``` r
# install.packages("devtools")
devtools::install_github("anthonydevaux/DynForest")
```

## Acknowledgements

We thank Dr. Louis Capitaine for FrechForest R code used in DynForest.
