SPAR
======================

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This R package enables you to apply Sparse Data-Driven Random Projection in Regression for High-Dimensional Data (see [Parzer, Filzmoser and Vana-Guer 2025](https://jdssv.org/index.php/jdssv/article/view/138)).
Exported functions are

- `spar`: performs the procedure for given thresholds lambda and numbers of marginal models, and acts as a help-function for the full cross-validated procedure spar.cv.
- `spar.cv`: performs the full procedure as described in the corresponding paper will cross-validation of the optimal threshold and number of models.

This repository includes development versions of CRAN package sparreg (https://github.com/lauravana/spareg), further development is moved to that package.

## Installation

```s
# install.packages("remotes")
remotes::install_github("RomanParzer/SPAR") # for current main branch
# remotes::install_github("RomanParzer/SPAR@*release") # for latest release
```

## Usage 

The method is designed for linear regression using a high-dimensional data set with more variables than observations.
The two main functions return objects, for which `coef`, `predict` and `plot` functions are available.

```s
data("example_data")
spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30))
spar_res
coefs <- coef(spar_res)
pred <- predict(spar_res,example_data$x)
plot(spar_res)
plot(spar_res,"MSE","nummod")
plot(spar_res,"numAct","lambda")
```

## License

This package is free and open source software, licensed under GPL-3.
