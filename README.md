
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
outcome in the paper: \> Devaux A., Helmer C., Dufouil C., Genuer R.,
Proust-Lima C. (2022). Random survival forests for competing risks with
multivariate longitudinal endogenous covariates. *arxiv*. [\<doi:
10.48550/arXiv.2208.05801\>](https://doi.org/10.48550/arXiv.2208.05801)

## Installation

You can install the development version of DynForest from
[GitHub](https://github.com/anthonydevaux/DynForest) with:

``` r
# install.packages("devtools")
devtools::install_github("anthonydevaux/DynForest")
#> Downloading GitHub repo anthonydevaux/DynForest@HEAD
#>          checking for file 'C:\Users\antho\AppData\Local\Temp\RtmpYDNKFY\remotes36e429807b74\anthonydevaux-DynForest-645e95b/DESCRIPTION' ...     checking for file 'C:\Users\antho\AppData\Local\Temp\RtmpYDNKFY\remotes36e429807b74\anthonydevaux-DynForest-645e95b/DESCRIPTION' ...   ✔  checking for file 'C:\Users\antho\AppData\Local\Temp\RtmpYDNKFY\remotes36e429807b74\anthonydevaux-DynForest-645e95b/DESCRIPTION' (342ms)
#>       ─  preparing 'DynForest': (440ms)
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'DynForest_0.1.0.tar.gz'
#>      
#> 
#> Installation du package dans 'C:/Users/antho/AppData/Local/Temp/RtmpglTtzB/temp_libpath413c4fa77ace'
#> (car 'lib' n'est pas spécifié)
```
