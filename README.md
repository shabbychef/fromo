

# fromo

[![Build Status](https://travis-ci.org/shabbychef/fromo.png)](https://travis-ci.org/shabbychef/fromo)
[![codecov.io](http://codecov.io/github/shabbychef/fromo/coverage.svg?branch=master)](http://codecov.io/github/shabbychef/fromo?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/fromo)](http://cran.rstudio.com/package=fromo) 
[![Downloads](http://cranlogs.r-pkg.org/badges/fromo?color=brightgreen)](http://www.r-pkg.org/pkg/fromo)

Fast robust moments via Rcpp, mostly as an exercise to learn Rcpp. 
Supports computation on vectors and matrices, and Monoidal append of moments (NYI).

-- Steven E. Pav, shabbychef@gmail.com

## Installation


```r
# get snapshot from github (may be buggy)
if (require(devtools)) {
    install_github("shabbychef/fromo")
}
```

# Basic Usage

## summary moments

Here is a speed comparison of the basic moment computation


```r
require(fromo)
require(moments)
require(microbenchmark)

x <- rnorm(1000)

dumbk <- function(x) {
    c(kurtosis(x) - 3, skewness(x), sd(x), mean(x), 
        length(x))
}

microbenchmark(kurt5(x), skew4(x), sd3(x), dumbk(x), 
    dumbk(x), kurtosis(x), skewness(x), sd(x), mean(x))
```

```
## Unit: microseconds
##         expr   min    lq  mean median    uq max neval
##     kurt5(x) 145.0 146.5 147.9  147.7 148.5 162   100
##     skew4(x)  83.9  85.3  87.1   86.0  87.0 120   100
##       sd3(x)  10.4  11.3  13.3   11.7  12.5 144   100
##     dumbk(x) 199.2 209.9 214.7  212.0 218.1 327   200
##  kurtosis(x)  86.9  91.6  95.0   92.6  93.9 256   100
##  skewness(x)  87.4  92.3  94.4   93.4  95.1 115   100
##        sd(x)  15.6  18.4  20.3   19.2  21.3  33   100
##      mean(x)   3.7   4.3   5.1    4.7   5.2  18   100
```

```r
x <- rnorm(1e+07, mean = 1e+12)

microbenchmark(kurt5(x), skew4(x), sd3(x), dumbk(x), 
    kurtosis(x), skewness(x), sd(x), mean(x), times = 10L)
```

```
## Unit: milliseconds
##         expr  min   lq mean median   uq  max neval
##     kurt5(x) 1449 1456 1496   1472 1498 1623    10
##     skew4(x)  820  821  830    826  831  854    10
##       sd3(x)   85   85   86     85   85   93    10
##     dumbk(x) 1718 1722 1769   1762 1808 1848    10
##  kurtosis(x)  841  843  869    858  883  930    10
##  skewness(x)  805  814  832    820  836  915    10
##        sd(x)   49   49   51     50   51   57    10
##      mean(x)   17   17   17     17   17   18    10
```

## running moments

Since an online algorithm is used, we can compute cumulative running moments. Moreover, we can 
//remove// observations, and thus compute moments over a fixed length lookback window. A demonstration:


```r
require(fromo)
require(moments)
require(microbenchmark)

set.seed(1234)
x <- rnorm(20)

k5 <- run_kurt5(x, winsize = 10L)
colnames(k5) <- c("excess_kurtosis", "skew", "stdev", 
    "mean", "nobs")
k5
```

```
##       excess_kurtosis  skew stdev   mean nobs
##  [1,]             NaN   NaN   NaN -1.207    1
##  [2,]             NaN   NaN  1.05 -0.465    2
##  [3,]             NaN -0.34  1.16  0.052    3
##  [4,]          -1.520 -0.13  1.53 -0.548    4
##  [5,]          -1.254 -0.50  1.39 -0.352    5
##  [6,]          -0.860 -0.79  1.30 -0.209    6
##  [7,]          -0.714 -0.70  1.19 -0.261    7
##  [8,]          -0.525 -0.64  1.11 -0.297    8
##  [9,]          -0.331 -0.58  1.04 -0.327    9
## [10,]          -0.331 -0.42  1.00 -0.383   10
## [11,]           0.262 -0.65  0.95 -0.310   10
## [12,]           0.017 -0.30  0.95 -0.438   10
## [13,]           0.699 -0.61  0.79 -0.624   10
## [14,]          -0.939  0.69  0.53 -0.383   10
## [15,]          -0.296  0.99  0.64 -0.330   10
## [16,]           1.078  1.33  0.57 -0.391   10
## [17,]           1.069  1.32  0.57 -0.385   10
## [18,]           0.868  1.29  0.60 -0.421   10
## [19,]           0.799  1.31  0.61 -0.449   10
## [20,]           1.193  1.50  1.07 -0.118   10
```

```r
# trust but verify
alt5 <- sapply(seq_along(x), function(iii) {
    rowi <- max(1, iii - 10 + 1)
    kurtosis(x[rowi:iii]) - 3
}, simplify = TRUE)

cbind(alt5, k5[, 1])
```

```
##         alt5       
##  [1,]    NaN    NaN
##  [2,] -2.000    NaN
##  [3,] -1.500    NaN
##  [4,] -1.520 -1.520
##  [5,] -1.254 -1.254
##  [6,] -0.860 -0.860
##  [7,] -0.714 -0.714
##  [8,] -0.525 -0.525
##  [9,] -0.331 -0.331
## [10,] -0.331 -0.331
## [11,]  0.262  0.262
## [12,]  0.017  0.017
## [13,]  0.699  0.699
## [14,] -0.939 -0.939
## [15,] -0.296 -0.296
## [16,]  1.078  1.078
## [17,]  1.069  1.069
## [18,]  0.868  0.868
## [19,]  0.799  0.799
## [20,]  1.193  1.193
```

