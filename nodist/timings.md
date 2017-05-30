

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
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.2 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] methods stats   utils   base   
## 
## other attached packages:
##  [1] RcppParallel_4.3.20    RcppRoll_0.2.2         roll_1.0.7             RollingWindow_0.2      microbenchmark_1.4-2.1
##  [6] moments_0.14           dplyr_0.5.0            fromo_0.1.3.3200       ggplot2_2.2.1          knitr_1.15.1          
## [11] devtools_1.12.0        Quandl_2.8.0           xts_0.9-7              zoo_1.7-12             drat_0.1.2            
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.10     magrittr_1.5     grDevices_3.3.2  munsell_0.4.3    colorspace_1.3-2 lattice_0.20-33  R6_2.2.0        
##  [8] stringr_1.2.0    httr_1.2.1       plyr_1.8.4       tools_3.3.2      grid_3.3.2       gtable_0.2.0     DBI_0.6-1       
## [15] withr_1.0.2      assertthat_0.2.0 lazyeval_0.2.0   digest_0.6.12    tibble_1.3.0     formatR_1.5      graphics_3.3.2  
## [22] memoise_1.1.0    evaluate_0.10    stringi_1.1.5    scales_0.4.1     jsonlite_1.4
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
##         expr   min    lq  mean median    uq   max neval   cld
##       sum(x)    80    80    85     81    87   151   100 a    
##      mean(x)   162   166   180    172   187   287   100 a    
##        sd(x)   459   469   506    487   513   806   100 ab   
##  skewness(x)  7872  8062  8548   8174  8643 12910   100   c  
##  kurtosis(x)  7832  7953  8392   8093  8620 10588   100   c  
##       sd3(x)   855   868   924    880   950  1248   100  b   
##     skew4(x)  7937  8121  8425   8223  8444 10818   100   c  
##     kurt5(x) 14262 14403 14989  14642 15291 21420   100    d 
##     dumbk(x) 16421 16839 18388  17670 18436 62493   100     e
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 15380 15612 15917  15724 16030 19572   100
##                                                               sd3(x, wts = w)   982   989  1047   1023  1075  1289   100
##                                                                 slow_sd(x, w)  1386  1472  3929   2056  2885 97461   100
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
##  as.centsums(x1, 1)  182  184  194    187  198  321   100  b  
##  as.centsums(x1, 2)  115  117  128    118  127  380   100 a   
##  as.centsums(x1, 3)  810  813  852    825  881 1208   100   c 
##  as.centsums(x1, 4) 1400 1405 1449   1420 1451 1795   100    d
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
##  c(obj1, obj2)  15 17   27     20 30 204   100   a
##  obj3 %-% obj1  12 14   24     17 28 121   100   a
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
##                        expr min  lq mean median  uq  max neval cld
##    join_cent_sums(rs1, rs2) 2.3 2.4  3.1    2.6 2.8 33.8   100   a
##  unjoin_cent_sums(rs3, rs2) 2.1 2.2  2.5    2.4 2.6  4.4   100   a
##  unjoin_cent_sums(rs3, rs1) 2.1 2.2  3.0    2.4 2.5 44.7   100   a
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
##  as.centcosums(x1, max_ord)  51 54   60     55 65 155   100   b
##             mobj3 %-% mobj1  17 18   21     19 20  39   100  a
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
        wins, robust = TRUE), running_sum(x, wins, 
        robust = FALSE), running_mean(x, wins), roll::roll_sum(xm, 
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
##                                                           expr   min    lq  mean median    uq    max neval      cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 33426 34551 39205  36082 37389 114737   100       g 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 66640 69059 73449  70528 72573 143229   100        h
##                            running_sum(x, wins, robust = TRUE)   146   151   160    154   167    205   100 a       
##                           running_sum(x, wins, robust = FALSE)    26    30    53     32    41   1676   100 a       
##                                          running_mean(x, wins)   147   151   184    158   169   2238   100 a       
##                                       roll::roll_sum(xm, wins)  1915  1964  2098   2032  2140   4231   100  b d    
##                                      roll::roll_mean(xm, wins)  2068  2112  2301   2206  2369   3531   100   cd    
##                                        roll::roll_sd(xm, wins)  5725  5842  6180   5926  6081   8346   100     e   
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   347   369   491    397   461   2579   100 abc     
##                             RollingWindow::RollingSum(x, wins)   108   128   259    155   195   3404   100 ab      
##                            RollingWindow::RollingMean(x, wins)   139   156   301    172   218   2643   100 ab      
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1737  1928  2044   1968  2099   4203   100  b d    
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1753  1978  2105   2033  2123   4119   100  b d    
##                                         running_skew4(x, wins)  3442  3547  3855   3615  4057   5805   100    d    
##                                         running_kurt5(x, wins)  5507  5648  6106   5735  6181  10824   100     e   
##                                         running_tstat(x, wins)   664   677   728    705   735   1038   100 abc     
##                                       running_zscored(x, wins)   664   676   729    696   751   1055   100 abc     
##                                        running_sharpe(x, wins)   668   680   775    702   775   3008   100 abc     
##                                    running_apx_median(x, wins) 13105 13362 14083  13563 14316  20137   100      f  
##                                      running_centered(x, wins)   550   559   606    578   619   1103   100 abc     
##                                        running_scaled(x, wins)   661   671   719    687   724   1057   100 abc
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

