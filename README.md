
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

> Devaux A., Helmer C., Genuer R. and Proust-Lima C. (2023). Random
> survival forests with multivariate longitudinal endogenous covariates.
> *Statistical Methods in Medical Research*.
> [\<doi:10.1177/09622802231206477\>](https://doi.org/10.1177/09622802231206477)

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

## Quick example

### Manage data

``` r
library(DynForest)
#> Registered S3 method overwritten by 'cmprsk':
#>   method      from
#>   plot.cuminc lcmm
data(pbc2)

# Get Gaussian distribution for longitudinal predictors
pbc2$serBilir <- log(pbc2$serBilir)
pbc2$SGOT <- log(pbc2$SGOT)
pbc2$albumin <- log(pbc2$albumin)
pbc2$alkaline <- log(pbc2$alkaline)
```

### Build DynForest objects

``` r
# Build longitudinal data
timeData <- pbc2[,c("id","time",
                    "serBilir","SGOT",
                    "albumin","alkaline")]

# Create object with longitudinal association for each predictor
timeVarModel <- list(serBilir = list(fixed = serBilir ~ time,
                                     random = ~ time),
                     SGOT = list(fixed = SGOT ~ time + I(time^2),
                                 random = ~ time + I(time^2)),
                     albumin = list(fixed = albumin ~ time,
                                    random = ~ time),
                     alkaline = list(fixed = alkaline ~ time,
                                     random = ~ time))
# Build fixed data
fixedData <- unique(pbc2[,c("id","age","drug","sex")])

# Build outcome data
Y <- list(type = "surv",
          Y = unique(pbc2[,c("id","years","event")]))
```

### Run `DynForest()` function

``` r
# Run DynForest function
res_dyn <- DynForest(timeData = timeData, fixedData = fixedData,
                     timeVar = "time", idVar = "id",
                     timeVarModel = timeVarModel, Y = Y,
                     ntree = 50, nodesize = 5, minsplit = 5,
                     cause = 2, ncores = 2, seed = 1234)
```

### Get summary

``` r
summary(res_dyn)
#> DynForest executed for survival (competing risk) outcome 
#>  Splitting rule: Fine & Gray statistic test 
#>  Out-of-bag error type: Integrated Brier Score 
#>  Leaf statistic: Cumulative incidence function 
#> ---------------- 
#> Input 
#>  Number of subjects: 312 
#>  Longitudinal: 4 predictor(s) 
#>  Numeric: 1 predictor(s) 
#>  Factor: 2 predictor(s) 
#> ---------------- 
#> Tuning parameters 
#>  mtry: 3 
#>  nodesize: 5 
#>  minsplit: 5 
#>  ntree: 50 
#> ---------------- 
#> ---------------- 
#> DynForest summary 
#>  Average depth per tree: 5.94 
#>  Average number of leaves per tree: 20.44 
#>  Average number of subjects per leaf: 9.67 
#>  Average number of events of interest per leaf: 4.34 
#> ---------------- 
#> Computation time 
#>  Number of cores used: 2 
#>  Time difference of 8.185554 mins
#> ----------------
```
