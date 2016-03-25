

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

For now, check this out:


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
##         expr   min    lq  mean median    uq max neval   cld
##     kurt5(x) 139.6 141.3 149.4  142.8 144.6 210   100    d 
##     skew4(x)  81.5  83.2  89.9   84.3  85.9 176   100   c  
##       sd3(x)  18.4  19.5  21.9   19.8  20.9  39   100  b   
##     dumbk(x) 194.4 197.5 218.8  200.3 211.9 324   200     e
##  kurtosis(x)  85.5  87.2  98.7   89.0  96.4 160   100   c  
##  skewness(x)  85.4  87.3  93.4   88.1  89.4 163   100   c  
##        sd(x)  14.5  16.2  19.1   17.1  18.2  52   100  b   
##      mean(x)   3.8   4.2   5.2    4.6   5.3  13   100 a
```

```r
x <- rnorm(1e+07, mean = 1e+12)

microbenchmark(kurt5(x), skew4(x), sd3(x), dumbk(x), 
    kurtosis(x), skewness(x), sd(x), mean(x), times = 10L)
```

```
## Unit: milliseconds
##         expr  min   lq mean median   uq  max neval    cld
##     kurt5(x) 1392 1394 1454   1415 1466 1687    10     e 
##     skew4(x)  801  808  817    814  832  836    10   c   
##       sd3(x)  168  169  170    170  172  175    10  b    
##     dumbk(x) 2332 2393 2509   2523 2594 2738    10      f
##  kurtosis(x) 1170 1174 1206   1192 1206 1330    10    d  
##  skewness(x) 1116 1142 1220   1230 1290 1353    10    d  
##        sd(x)   46   47   49     48   49   55    10 a     
##      mean(x)   16   16   16     16   17   17    10 a
```

