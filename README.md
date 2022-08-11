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



# TRexSelector
**Title**: The T-Rex selector for fast high-dimensional variable selection with FDR control

**Description**: It performs fast variable selection in large-scale high-dimensional settings while controlling the false discovery rate (FDR) at a user-defined target level. The package is based on the T-Rex selector paper (available at https://arxiv.org/abs/2110.06048).

**Note**: The T-Rex selector performs terminated-random experiments (T-Rex) using the T-LARS algorithm ([R package](https://CRAN.R-project.org/package=tlars)) and fuses the selected active sets of all random experiments to obtain a final set of selected variables. The T-Rex selector provably controls the false discovery rate (FDR), i.e., the expected fraction of selected false positives among all selected variables, at the user-defined target level while maximizing the number of selected variables and, thereby, achieving a high true positive rate (TPR) (i.e., power). The T-Rex selector can be applied in various fields, such as genomics, financial engineering, or any other field that requires a fast and FDR-controlling variable/feature selection method for large-scale high-dimensional settings.

In the following sections, we show you how to install and use the package.

## Installation
Before installing the 'TRexSelector' package, you need to install the required 'tlars' package. You can install the 'tlars' package from [CRAN](https://CRAN.R-project.org/package=tlars) or [GitHub](https://github.com/jasinmachkour/tlars) with: 


```r
# Install stable version from CRAN
install.packages("tlars")

# Install development version from GitHub
install.packages("devtools")
devtools::install_github("jasinmachkour/tlars")
```

Then, you can install the 'TRexSelector' package with:


```r
devtools::install_github("jasinmachkour/TRexSelector")
```

You can open the help pages with:

```r
library(TRexSelector)
help(package = "TRexSelector")
?trex
?random_experiments
?lm_dummy
?add_dummies
?add_dummies_GVS
?FDP
?TPP
# etc.
```

To cite the package 'TRexSelector' in publications use:

```r
citation("TRexSelector")
```

# Quick Start
This section illustrates the basic usage of the 'TRexSelector' package to perform FDR-controlled variable selection in large-scale high-dimensional settings based on the T-Rex selector.

1. **First**, we generate a high-dimensional Gaussian data set with sparse support:


```r
library(TRexSelector)

# Setup
n <- 75 # number of observations
p <- 150 # number of variables
num_act <- 3 # number of true active variables
beta <- c(rep(1, times = num_act), rep(0, times = p - num_act)) # coefficient vector
true_actives <- which(beta > 0) # indices of true active variables
num_dummies <- p # number of dummy predictors (also referred to as dummies)

# Generate Gaussian data
set.seed(123)
X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y <- X %*% beta + stats::rnorm(n)
```

2. **Second**, we perform FDR-controlled variable selection using the T-Rex selector for a target FDR of 5%:


```r
# Seed
set.seed(1234)

# Numerical zero
eps <- .Machine$double.eps

# Variable selection via T-Rex
res <- trex(X = X, y = y, tFDR = 0.05, verbose = FALSE)
selected_var <- which(res$selected_var > eps)
paste0("True active variables: ", paste(as.character(true_actives), collapse = ", "))
#> [1] "True active variables: 1, 2, 3"
paste0("Selected variables: ", paste(as.character(selected_var), collapse = ", "))
#> [1] "Selected variables: 1, 2, 3"
```

So, for a preset target FDR of 5%, the T-Rex selector has selected all true active variables and there are no false positives in this example.

Note that users have to choose the target FDR according to the requirements of their specific applications.

## Documentation
For more information and some examples, please check the [GitHub-vignette](https://htmlpreview.github.io/?https://github.com/jasinmachkour/TRexSelector/blob/main/vignettes/trex_usage_and_simulations.html).

## Links
T-Rex paper: https://arxiv.org/abs/2110.06048

'TRexSelector' package: [GitHub-TRexSelector](https://github.com/jasinmachkour/TRexSelector).

README file: [GitHub-readme](https://htmlpreview.github.io/?https://github.com/jasinmachkour/TRexSelector/blob/main/README.html).

Vignette: [GitHub-vignette](https://htmlpreview.github.io/?https://github.com/jasinmachkour/TRexSelector/blob/main/vignettes/trex_usage_and_simulations.html).

'tlars' package: [CRAN-tlars](https://CRAN.R-project.org/package=tlars) and [GitHub-tlars](https://github.com/jasinmachkour/tlars).


