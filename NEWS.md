# DynForest 1.1.3 (CRAN version)

-   Improved `DynForest()` computation time function with survival outcome

# DynForest 1.1.2 (CRAN version)

-   Added `conv_issue` return object for each tree providing information regarding convergence of the longitudinal data
-   Added S3 method `print()` for class `DynForest`, `DynForestVIMP`, `DynForestgVIMP` and `DynForestVarDepth`
-   Added S3 method `plot()` to replace `plot_nodeCIF()` and `plot_CIF()`

# DynForest 1.1.1

-   Fixed bugs with `Inputs` object in `compute_VIMP()` and `compute_gVIMP()` functions (thanks to ChenXinaha GitHub user)
-   Fixed an issue with single event in survival mode (thanks to Dr. Labruère Chazal)
-   Fixed bug with time variable name in formula (thanks to Dr. Labruère Chazal)
-   Added `getTree()` function to extract split information from a tree specified by user
-   Added `getTreeNodes()` function to extract nodes identifiers for a given tree
-   Added `plot_nodeCIF()` function to plot estimated CIF for given tree nodes

# DynForest 1.1.0

-   Added S3 method `plot()` to replace `plot_mindepth()`, `plot_VIMP()` and `plot_gVIMP()` functions
-   Added a vignette to introduce `DynForest` methodology
-   Added a vignette for each outcome (survival, categorical and continuous)
-   Output names `Curve` and `Scalar` were replaced by `Longitudinal` and `Numeric`, respectively
-   Improved permutation process in `compute_VIMP()` and `compute_gVIMP()` functions
