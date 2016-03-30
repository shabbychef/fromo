

# fromo

[![Build Status](https://travis-ci.org/shabbychef/fromo.png)](https://travis-ci.org/shabbychef/fromo)
[![codecov.io](http://codecov.io/github/shabbychef/fromo/coverage.svg?branch=master)](http://codecov.io/github/shabbychef/fromo?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/fromo)](http://cran.rstudio.com/package=fromo) 
[![Downloads](http://cranlogs.r-pkg.org/badges/fromo?color=brightgreen)](http://www.r-pkg.org/pkg/fromo)

Fast robust moments via Rcpp, mostly as an exercise to learn Rcpp. 
Supports computation on vectors and matrices, and Monoidal append (and unappend) of moments.

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
##         expr   min    lq mean median    uq    max neval cld
##     kurt5(x) 184.6 199.0  223  214.5 247.4    304   100   a
##     skew4(x) 109.0 121.3  136  130.8 146.8    291   100   a
##       sd3(x)  14.5  17.2   20   19.1  21.5     43   100   a
##     dumbk(x) 279.8 305.7  940  337.1 371.3 120142   200   a
##  kurtosis(x) 121.5 130.7  146  145.0 158.7    193   100   a
##  skewness(x) 121.8 132.4  149  146.1 164.2    194   100   a
##        sd(x)  24.6  28.9   33   31.7  36.6     57   100   a
##      mean(x)   5.2   6.8    8    7.6   8.5     22   100   a
```

```r
x <- rnorm(1e+07, mean = 1e+12)

microbenchmark(kurt5(x), skew4(x), sd3(x), dumbk(x), 
    kurtosis(x), skewness(x), sd(x), mean(x), times = 10L)
```

```
## Unit: milliseconds
##         expr  min   lq mean median   uq  max neval     cld
##     kurt5(x) 1815 1817 1820   1818 1821 1830    10      f 
##     skew4(x) 1047 1047 1054   1054 1059 1067    10    d   
##       sd3(x)  101  101  101    101  101  103    10   c    
##     dumbk(x) 2194 2197 2205   2201 2208 2239    10       g
##  kurtosis(x) 1083 1083 1091   1085 1086 1135    10     e  
##  skewness(x) 1035 1036 1044   1037 1038 1108    10    d   
##        sd(x)   58   58   58     58   58   59    10  b     
##      mean(x)   19   19   19     19   19   20    10 a
```

## running moments

Since an online algorithm is used, we can compute cumulative running moments. Moreover, we can 
//remove// observations, and thus compute moments over a fixed length lookback window. The code
checks for negative even moments caused by roundoff, and restarts the computation to correct;
periodic recomputation can be forced by an input parameter.

A demonstration:


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

## running 'scale'

Through template magic, the same code can perform running centering, scaling, z-scoring and so on:


```r
require(fromo)
require(moments)
require(microbenchmark)

set.seed(1234)
x <- rnorm(20)

xz <- run_zscored(x, winsize = 10L)

# trust but verify
altz <- sapply(seq_along(x), function(iii) {
    rowi <- max(1, iii - 10 + 1)
    (x[iii] - mean(x[rowi:iii]))/sd(x[rowi:iii])
}, simplify = TRUE)

cbind(xz, altz)
```

```
##              altz
##  [1,]   NaN    NA
##  [2,]  0.71  0.71
##  [3,]  0.89  0.89
##  [4,] -1.18 -1.18
##  [5,]  0.56  0.56
##  [6,]  0.55  0.55
##  [7,] -0.26 -0.26
##  [8,] -0.23 -0.23
##  [9,] -0.23 -0.23
## [10,] -0.51 -0.51
## [11,] -0.17 -0.17
## [12,] -0.59 -0.59
## [13,] -0.19 -0.19
## [14,]  0.84  0.84
## [15,]  2.02  2.02
## [16,]  0.49  0.49
## [17,] -0.22 -0.22
## [18,] -0.82 -0.82
## [19,] -0.64 -0.64
## [20,]  2.37  2.37
```
