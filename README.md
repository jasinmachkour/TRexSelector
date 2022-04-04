
# tknock

**Title**: The T-Knock filter for fast high-dimensional variable
selection with FDR control

**Description**: It performs variable selection in high-dimensional
settings where the number of variables is potentially higher than the
number of observations (data points) while controlling the FDR at a
user-defined target level. The T-Knock filter paper is available at
<https://arxiv.org/abs/2110.06048>.

## Installation

Before installing the tknock package, you need to install the required
tlars package with:

``` r
install.packages("devtools")
devtools::install_github("jasinmachkour/tlars")
```

You can install the tknock package with:

``` r
devtools::install_github("jasinmachkour/tknock")
```
