
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/anthonydevaux/DynForest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anthonydevaux/DynForest/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Description

`DynForest` is a R package aiming to predict an outcome using
multivariate longitudinal predictors. The method is based on random
forest principle where the longitudinal predictors are modeled through
the random forest. `DynForest` currently supports continuous, discrete
and survival outcome. The methodology is fully described for a survival
outcome in the paper:

> Devaux A., Helmer C., Dufouil C., Genuer R., Proust-Lima C. (2022).
> Random survival forests for competing risks with multivariate
> longitudinal endogenous covariates. *arXiv*. [\<doi:
> 10.48550/arXiv.2208.05801\>](https://doi.org/10.48550/arXiv.2208.05801)

## Installation

You can install the development version of DynForest from
[GitHub](https://github.com/anthonydevaux/DynForest) with:

``` r
# install.packages("devtools")
devtools::install_github("anthonydevaux/DynForest")
```
