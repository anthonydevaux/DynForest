# DynForest 1.1.1

* Fixed bugs with `Inputs` object in `compute_VIMP()` and `compute_gVIMP()` functions (thanks to ChenXinaha GitHub user)
* Fixed an issue with single event in survival mode (thanks to Dr. Labruère Chazal)

# DynForest 1.1.0

* Added S3 method `plot()` to replace `plot_mindepth()`, `plot_VIMP()` and `plot_gVIMP()` functions
* Added a vignette to introduce `DynForest` methodology 
* Added a vignette for each outcome (survival, categorical and continuous)
* Output names `Curve` and `Scalar` were replaced by `Longitudinal` and `Numeric`, respectively
* Improved permutation process in `compute_VIMP()` and `compute_gVIMP()` functions
