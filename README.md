
# HIDeM

<!-- badges: start -->
<!-- badges: end -->

The goal of HIDeM is to perform variable selection on illness-death models with interval-censored data. 

## Installation

You can install the development version of HIDeM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library("devtools")
install_github("arianebercu/HIDeM",build_vignettes=T)
```

## Examples

In this package, we choose to illustrate a simulated scenario and an application to Paquid data set using vignettes: Reg-IDM-ICT-example and Reg-IDM-ICT-Paq1000. The Reg-IDM-ICT-example takes back the simulated Scenario B -see paper- and presents the variable selection process. 

👉 See the [Getting Started vignette](doc/Reg-IDM-ICT-example.html and doc/Reg-IDM-ICT-Paq1000.html)

``` r
library("HIDeM")
vignette("Reg-IDM-ICT-example", package = "HIDeM")
vignette("Reg-IDM-ICT-Paq1000", package = "HIDeM")
## basic example code
```

