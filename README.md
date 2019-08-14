## README

This is an ad hoc package of quantitative genetics functions implemented using the C++ [Eigen matrix algebra library](<http://eigen.tuxfamily.org/index.php?title=Main_Page>). These are intended to speed up computations using high-throughput genotyping data. Many of these functions can be found in other software or R packages. Re-implementing them here either allows for increased flexibility in processing output or increasing computational efficiency.

The functions assume that marker data are stored in R matrices. I currently have no plans to provide support for out-of-memory matrices at this time.

These functions are written and updated as necessary for my own research purposes.

### Installation

Please note that this package depends on the [FastMath](<https://github.com/amkusmec/FastMath>) package for reimplementations of commonly used mathematical functions.

```
install.packages(c("devtools", "Rcpp", "RcppEigen", "mgcv"))
devtools::install_github("amkusmec/FastMath")
devtools::install_github("amkusmec/QGenTools")
note: Rcpp v.1.0.2 and higher is required to finish installation of FastMath. Currently CRAN only provide earlier version.
```

### Performance

These functions are not strongly optimized, mostly relying on the speed advantages of C++ over R and internal optimization by the Eigen library during compilation. You may derive greater performance by manipulating compiler flags, but the R defaults should perform well enough for most situations. On the other hand, these functions are not optimized at all with respect to memory usage. Care should be taken when requesting linkage disequilibrium, for example, on large marker matrices because all pairwise correlations will be returned.
