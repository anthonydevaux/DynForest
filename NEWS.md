# DynForest 1.1.1

* Fixed bugs with `Inputs` object in `compute_VIMP()` and `compute_gVIMP()` functions (thanks to ChenXinaha GitHub user)
* Fixed an issue with single event in survival mode (thanks to Dr. Labruère Chazal)
* Fixed bug with time variable name in formula (thanks to Dr. Labruère Chazal)
* Added `getTree()` function to extract split information from a tree specified by user
* Added `plot_nodeCIF()` function to plot estimated CIF for given tree nodes

# DynForest 1.1.0

* Added S3 method `plot()` to replace `plot_mindepth()`, `plot_VIMP()` and `plot_gVIMP()` functions
* Added a vignette to introduce `DynForest` methodology 
* Added a vignette for each outcome (survival, categorical and continuous)
* Output names `Curve` and `Scalar` were replaced by `Longitudinal` and `Numeric`, respectively
* Improved permutation process in `compute_VIMP()` and `compute_gVIMP()` functions
