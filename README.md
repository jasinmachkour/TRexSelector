---
output:
  html_document:
    variant: markdown_github
    keep_md: yes
  md_document:
    variant: markdown_github
  pdf_document: default
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# tknock
**Title**: The T-Knock filter for fast high-dimensional variable selection with FDR control

**Description**: It performs fast variable selection in large-scale high-dimensional settings while controlling the false discovery rate (FDR) at a user-defined target level. The package is based on the T-Knock filter paper (available at https://arxiv.org/abs/2110.06048).

**Note**: The T-Rex selector performs terminated-random experiments (T-Rex) using the T-LARS algorithm ([R package](https://github.com/jasinmachkour/tlars)) and fuses the selected active sets of all random experiments to obtain a final set of selected variables. The T-Rex selector provably controls the false discovery rate (FDR), i.e., the expected fraction of selected false positives among all selected variables, at the user-defined target level while maximizing the number of selected variables and, thereby, achieving a high true positive rate (TPR) (i.e., power). The T-Rex selector can be applied in various fields, such as genomics, financial engineering, or any other field that requires a fast and FDR-controlling variable/feature selection method for large-scale high-dimensional settings.

In the following, we show how to use the package and show you how genome-wide association studies (GWAS) can be performed using the T-Rex selector. 

## Installation

Before installing the tknock package, you need to install the required tlars package. You can install the tlars package from [GitHub](https://github.com/jasinmachkour/tlars) with 

``` r
install.packages("devtools")
devtools::install_github("jasinmachkour/tlars")
```

Then, you can install the tknock package with:

``` r
devtools::install_github("jasinmachkour/tknock")
```

You can open the help pages with

```r
library(tknock)
help(package = "tknock")
?tknock
?random_experiments
?lm_dummy
?add_dummies
?add_dummies_GVS
?FDP
?TPP
# etc.
```

To cite the package ‘tknock’ in publications use:

```r
citation("tknock")
```

# Quick Start
In the following, we illustrate the basic usage of the tknock package to perform FDR-controlled variable selection in large-scale high-dimensional settings using the T-Rex selector:

1. **First**, we generate a high-dimensional Gaussian data set with sparse support:


```r
library(tknock)

# Setup
n <- 75 # Number of observations
p <- 150 # Number of variables
num_act <- 3 # Number of true active variables
beta <- c(rep(1, times = num_act), rep(0, times = p - num_act)) # Coefficient vector
true_actives <- which(beta > 0) # Indices of true active variables
num_dummies <- p # Number of dummy predictors (or dummies)

# Generate Gaussian data
set.seed(123)
X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y <- X %*% beta + stats::rnorm(n)
```

2. **Second**, we perform FDR-controlled variable selection using the T-Rex selector with a target FDR of 5%:


```r
# Seed
set.seed(1234)

# Numerical zero
eps <- .Machine$double.eps

# Variable selection via T-Rex
res <- tknock(X = X, y = y, tFDR = 0.05, verbose = FALSE)
selected_var <- which(res$selected_var > eps)
paste0("True active variables: ", paste(as.character(true_actives), collapse = ", "))
#> [1] "True active variables: 1, 2, 3"
paste0("Selected variables: ", paste(as.character(selected_var), collapse = ", "))
#> [1] "Selected variables: 1, 2, 3"
```

So, with a target FDR of 5%, the T-Rex selector has selected all true active variables and there is no false positive in this example.

Note that users have to specify the target FDR according to the requirements of their specific applications.

## Documentation
For more information and some examples, please check the ...

## Links
T-Rex paper: https://arxiv.org/abs/2110.06048

trex package: [GitHub-trex](https://github.com/jasinmachkour/tknock).

README file: 

Vignette:

tlars package: [GitHub-tlars](https://github.com/jasinmachkour/tlars).
