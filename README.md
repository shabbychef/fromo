

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
##         expr   min    lq  mean median    uq    max neval  cld
##     kurt5(x) 186.8 189.4 197.1  191.1 196.5  413.8   100   c 
##     skew4(x) 111.3 113.5 120.6  115.2 116.9  403.6   100  b  
##       sd3(x)  13.7  15.2  17.1   16.3  17.1   69.1   100 a   
##     dumbk(x) 265.2 270.5 294.5  279.3 289.7 1543.4   200    d
##  kurtosis(x) 115.2 118.1 134.5  120.7 125.4 1233.1   100  b  
##  skewness(x) 116.3 118.6 128.8  122.4 125.3  398.0   100  b  
##        sd(x)  22.4  24.8  27.3   26.7  28.0   45.3   100 a   
##      mean(x)   5.6   6.2   6.9    6.8   7.6    8.9   100 a
```

```r
x <- rnorm(1e+07, mean = 1e+12)

microbenchmark(kurt5(x), skew4(x), sd3(x), dumbk(x), 
    kurtosis(x), skewness(x), sd(x), mean(x), times = 10L)
```

```
## Unit: milliseconds
##         expr  min   lq mean median   uq  max neval      cld
##     kurt5(x) 1908 1910 1912   1911 1912 1926    10       g 
##     skew4(x) 1104 1105 1109   1105 1113 1126    10     e   
##       sd3(x)  104  104  105    105  105  109    10   c     
##     dumbk(x) 2271 2275 2284   2277 2295 2310    10        h
##  kurtosis(x) 1119 1121 1134   1124 1154 1159    10      f  
##  skewness(x) 1072 1072 1078   1074 1078 1113    10    d    
##        sd(x)   59   60   60     60   61   61    10  b      
##      mean(x)   20   20   20     20   20   21    10 a
```

