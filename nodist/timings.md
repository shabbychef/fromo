

# fromo timings

To prevent performance regressions, compare them here, including benchmarks
against other packages:



```r
library(fromo)
library(RollingWindow)
library(roll)
library(RcppRoll)
library(microbenchmark)
library(moments)
library(RcppParallel)
# keep this constant for comparison
setThreadOptions(numThreads = 2)

print(sessionInfo())
```

```
## R version 3.3.0 RC (2016-05-01 r70566)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux stretch/sid
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] utils   methods base   
## 
## other attached packages:
##  [1] RcppParallel_4.3.20    RcppRoll_0.2.2         roll_1.0.7             RollingWindow_0.2      microbenchmark_1.4-2.1
##  [6] moments_0.14           dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.3400      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11     grDevices_3.3.0  assertthat_0.1   R6_2.1.2         grid_3.3.0       plyr_1.8.4       DBI_0.6-1       
##  [8] gtable_0.2.0     formatR_1.4      magrittr_1.5     evaluate_0.10    scales_0.4.1     stringi_1.0-1    lazyeval_0.2.0  
## [15] graphics_3.3.0   tools_3.3.0      stringr_1.0.0    munsell_0.4.3    stats_3.3.0      colorspace_1.3-2 tibble_1.2
```

```r
set.seed(12345)
x <- rnorm(1e+05)
xpos <- runif(length(x)) + 1
xm <- matrix(x, ncol = 1)
xmps <- matrix(xpos, ncol = 1)
w <- runif(length(x))

dumbk <- function(x) {
    c(kurtosis(x) - 3, skewness(x), sd(x), mean(x), 
        length(x))
}

checkit <- microbenchmark(sum(x), mean(x), sd(x), skewness(x), 
    kurtosis(x), sd3(x), skew4(x), kurt5(x), dumbk(x))
print(checkit)
```

```
## Unit: microseconds
##         expr   min    lq  mean median    uq   max neval     cld
##       sum(x)    80    80    82     80    81   108   100 a      
##      mean(x)   162   164   168    165   168   283   100 a      
##        sd(x)   459   466   484    479   489   633   100  b     
##  skewness(x)  8306  8407  8662   8445  8542 10721   100     e  
##  kurtosis(x)  8126  8245  8461   8286  8373 10310   100    d   
##       sd3(x)   853   860   881    870   891  1126   100   c    
##     skew4(x)  8314  8408  8510   8462  8538  9252   100    de  
##     kurt5(x) 15053 15177 15315  15254 15366 17219   100      f 
##     dumbk(x) 17149 17321 17841  17451 18125 19607   100       g
```

```r
resdf <- checkit
```


```r
# weights
slow_sd <- function(x, w) {
    n0 <- length(x)
    mu <- weighted.mean(x, w = w)
    sg <- sqrt(sum(w * (x - mu)^2)/(n0 - 1))
    c(sg, mu, n0)
}
checkit <- microbenchmark(cent_moments(x, max_order = 4, 
    wts = w, na_rm = TRUE, normalize_wts = FALSE), 
    sd3(x, wts = w), slow_sd(x, w))
print(checkit)
```

```
## Unit: microseconds
##                                                                          expr   min    lq  mean median    uq   max neval
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16217 16339 16482  16409 16525 17892   100
##                                                               sd3(x, wts = w)   981   985  1014   1000  1026  1155   100
##                                                                 slow_sd(x, w)  1400  1415  1925   1481  2751  3802   100
##  cld
##    c
##  a  
##   b
```

```r
resdf <- rbind(resdf, checkit)
```


```r
set.seed(12345)
x1 <- runif(10000)
x2 <- runif(length(x1))

checkit <- microbenchmark(as.centsums(x1, 1), as.centsums(x1, 
    2), as.centsums(x1, 3), as.centsums(x1, 4))
print(checkit)
```

```
## Unit: microseconds
##                expr  min   lq mean median   uq  max neval  cld
##  as.centsums(x1, 1)  188  189  194    190  192  298   100  b  
##  as.centsums(x1, 2)  113  115  121    116  122  274   100 a   
##  as.centsums(x1, 3)  874  876  895    879  898 1040   100   c 
##  as.centsums(x1, 4) 1504 1510 1539   1523 1545 1671   100    d
```

```r
resdf <- rbind(resdf, checkit)

# join them together
max_ord <- 6L
obj1 <- as.centsums(x1, max_ord)
obj2 <- as.centsums(x2, max_ord)
obj3 <- as.centsums(c(x1, x2), max_ord)

checkit <- microbenchmark(c(obj1, obj2), obj3 %-% obj1)
print(checkit)
```

```
## Unit: microseconds
##           expr min lq mean median uq max neval cld
##  c(obj1, obj2)  15 16   19     16 17 185   100   b
##  obj3 %-% obj1  12 13   15     13 14  99   100  a
```

```r
resdf <- rbind(resdf, checkit)

max_ord <- 6L
rs1 <- cent_sums(x1, max_ord)
rs2 <- cent_sums(x2, max_ord)
rs3 <- cent_sums(c(x1, x2), max_ord)

checkit <- microbenchmark(join_cent_sums(rs1, rs2), 
    unjoin_cent_sums(rs3, rs2), unjoin_cent_sums(rs3, 
        rs1))
print(checkit)
```

```
## Unit: microseconds
##                        expr min  lq mean median  uq max neval cld
##    join_cent_sums(rs1, rs2) 2.2 2.3  2.9    2.4 2.6  27   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.1  2.4    2.2 2.4  14   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.1  2.4    2.2 2.4  15   100   a
```

```r
resdf <- rbind(resdf, checkit)
```


```r
set.seed(54321)
x1 <- matrix(rnorm(100 * 4), ncol = 4)
x2 <- matrix(rnorm(100 * 4), ncol = 4)

max_ord <- 2L

# join them together
mobj1 <- as.centcosums(x1, max_ord)
mobj2 <- as.centcosums(x2, max_ord)
mobj3 <- as.centcosums(rbind(x1, x2), max_ord)
alt3 <- c(mobj1, mobj2)
# unjoin them, with this one weird operator:
alt2 <- mobj3 %-% mobj1
alt1 <- mobj3 %-% mobj2

checkit <- microbenchmark(as.centcosums(x1, max_ord), 
    mobj3 %-% mobj1)
print(checkit)
```

```
## Unit: microseconds
##                        expr min lq mean median uq max neval cld
##  as.centcosums(x1, max_ord)  51 53   59     54 61 145   100   b
##             mobj3 %-% mobj1  16 18   21     19 20  55   100  a
```

```r
resdf <- rbind(resdf, checkit)
```


```r
set.seed(4422)
x <- rnorm(10000)
xpos <- runif(length(x)) + 1
xm <- matrix(x, ncol = 1)
xmps <- matrix(xpos, ncol = 1)
w <- runif(length(x))

dumb_zscore <- function(x, window) {
    altz <- sapply(seq_along(x), function(iii) {
        rowi <- max(1, iii - window + 1)
        xrang <- x[rowi:iii]
        (x[iii] - mean(xrang))/sd(xrang)
    }, simplify = TRUE)
}

# run fun on each wins sized window...
silly_fun <- function(x, wins, fun, ...) {
    xout <- rep(NA, length(x))
    for (iii in seq_along(x)) {
        xout[iii] <- fun(x[max(1, iii - wins + 1):iii], 
            ...)
    }
    xout
}

wins <- 250

checkit <- microbenchmark(silly_fun(x, wins, sum, na.rm = FALSE), 
    silly_fun(x, wins, mean, na.rm = FALSE), running_sum(x, 
        wins), running_mean(x, wins), roll::roll_sum(xm, 
        wins), roll::roll_mean(xm, wins), roll::roll_sd(xm, 
        wins), RollingWindow::RollingSum(x, wins, na_method = "ignore"), 
    RollingWindow::RollingSum(x, wins), RollingWindow::RollingMean(x, 
        wins), RcppRoll::roll_sum(xm, n = wins, align = "right", 
        fill = NA), RcppRoll::roll_mean(xm, n = wins, 
        align = "right", fill = NA), running_skew4(x, 
        wins), running_kurt5(x, wins), running_tstat(x, 
        wins), running_zscored(x, wins), running_sharpe(x, 
        wins), running_apx_median(x, wins), running_centered(x, 
        wins), running_scaled(x, wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr   min    lq  mean median    uq    max neval       cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32577 32985 35022  34775 35277 109846   100        h 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 65264 66399 69919  67832 68898 143829   100         i
##                                           running_sum(x, wins)    73    77    85     80    89    155   100 a        
##                                          running_mean(x, wins)    73    77   109     80    89   2243   100 a        
##                                       roll::roll_sum(xm, wins)  1792  1818  1898   1877  1941   2333   100  bcd     
##                                      roll::roll_mean(xm, wins)  1957  1983  2058   2034  2097   2450   100    de    
##                                        roll::roll_sd(xm, wins)  5500  5604  5717   5667  5737   6547   100      f   
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   359   374   525    388   436   2670   100 a  d     
##                             RollingWindow::RollingSum(x, wins)   108   121   255    135   187   2524   100 ab       
##                            RollingWindow::RollingMean(x, wins)   136   154   291    169   215   2461   100 a c      
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1737  1830  1965   1946  1991   4063   100   cd     
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1720  1900  2000   1999  2042   4070   100    d     
##                                         running_skew4(x, wins)  3609  3642  3745   3684  3774   6114   100     e    
##                                         running_kurt5(x, wins)  5838  5915  6072   5975  6061   8318   100      f   
##                                         running_tstat(x, wins)   652   664   682    668   688    878   100 a  d     
##                                       running_zscored(x, wins)   661   667   737    674   695   3074   100 a  d     
##                                        running_sharpe(x, wins)   658   664   686    673   692    858   100 a  d     
##                                    running_apx_median(x, wins) 14017 14096 14402  14205 14429  16770   100       g  
##                                      running_centered(x, wins)   572   578   605    588   605    831   100 a  d     
##                                        running_scaled(x, wins)   656   665   706    672   692   2811   100 a  d
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
```

```
## Error in library(readr): there is no package called 'readr'
```

```r
readr::write_csv(resdf, "timings.csv")
```

```
## Error in loadNamespace(name): there is no package called 'readr'
```

