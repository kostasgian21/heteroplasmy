# heteroplasmy

`heteroplasmy` is an R package to explore different fits to the Kimura distribution using mtDNA heteroplasmy data, and estimate the standard error of the variance (and other quantities related to uncertainty estimation). This repository is a fork containing the source code for this package, although its original home is https://github.com/kostasgian21/heteroplasmy . This repo also contains R scripts that use the package to produce plots illustrating some issues with heteroplasmy fitting and analysis. These scripts are `kimura-issues-package.R` and `kimura-plots-package.R`.

To install the `heteroplasmy` package you will need the `devtools` R library. If you have this, you can ignore the first line below.

```
install.packages("devtools")
library("devtools")
devtools::install_github("kostasgian21/heteroplasmy")`
```

You can then load the `heteroplasmy` package with

`library("heteroplasmy")`

The documentation is currently incomplete but we are working on it; please contact us if you have any questions!

