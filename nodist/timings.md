

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

print(Sys.info())
```

```
##                                       sysname                                       release                                       version 
##                                       "Linux"                            "4.4.0-77-generic" "#98-Ubuntu SMP Wed Apr 26 08:34:02 UTC 2017" 
##                                      nodename                                       machine                                         login 
##                                "9320074115d9"                                      "x86_64"                                     "unknown" 
##                                          user                                effective_user 
##                                        "spav"                                        "spav"
```

```r
print(sessionInfo())
```

```
## R version 3.3.0 RC (2016-05-01 r70566)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux stretch/sid
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
##  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] utils   methods base   
## 
## other attached packages:
##  [1] RcppParallel_4.3.20    RcppRoll_0.2.2         roll_1.0.7             RollingWindow_0.2      microbenchmark_1.4-2.1 moments_0.14          
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4000      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11     grDevices_3.3.0  assertthat_0.1   R6_2.1.2         grid_3.3.0       plyr_1.8.4       DBI_0.6-1        gtable_0.2.0     formatR_1.4     
## [10] magrittr_1.5     evaluate_0.10    scales_0.4.1     stringi_1.0-1    lazyeval_0.2.0   graphics_3.3.0   tools_3.3.0      stringr_1.0.0    munsell_0.4.3   
## [19] stats_3.3.0      colorspace_1.3-2 tibble_1.2
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
##       sum(x)    80    80    83     80    81   143   100 a      
##      mean(x)   162   164   172    165   177   278   100 a      
##        sd(x)   459   463   483    471   490   720   100  b     
##  skewness(x)  8287  8362  8750   8443  8888 13124   100     e  
##  kurtosis(x)  8145  8200  8508   8267  8490 10415   100    d   
##       sd3(x)   852   859   883    868   898  1061   100   c    
##     skew4(x)  8307  8393  8586   8444  8629 10349   100    de  
##     kurt5(x) 14982 15078 15467  15203 15638 18416   100      f 
##     dumbk(x) 17063 17232 17953  17403 18684 21700   100       g
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
##                                                                          expr   min    lq  mean median    uq   max neval cld
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16104 16301 16915  16443 17643 21456   100   c
##                                                               sd3(x, wts = w)   981   985  1022   1000  1062  1133   100 a  
##                                                                 slow_sd(x, w)  1399  1422  1960   1512  2778  3986   100  b
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
##  as.centsums(x1, 1)  188  191  204    202  208  287   100  b  
##  as.centsums(x1, 2)  113  115  124    117  125  318   100 a   
##  as.centsums(x1, 3)  847  850  907    872  942 1149   100   c 
##  as.centsums(x1, 4) 1472 1477 1544   1502 1599 1935   100    d
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
##  c(obj1, obj2)  14 15   19     15 16 202   100   b
##  obj3 %-% obj1  11 11   13     12 13  93   100  a
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
##    join_cent_sums(rs1, rs2) 2.2 2.4  2.9    2.5 2.6  26   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.1  2.5    2.2 2.5  21   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.1  2.5    2.3 2.4  15   100   a
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
##  as.centcosums(x1, max_ord)  51 53   58     54 64 144   100   b
##             mobj3 %-% mobj1  16 17   20     18 20  32   100  a
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
        wins), RollingWindow::RollingStd(x, wins), 
    RcppRoll::roll_sum(xm, n = wins, align = "right", 
        fill = NA), RcppRoll::roll_mean(xm, n = wins, 
        align = "right", fill = NA), RcppRoll::roll_sd(xm, 
        n = wins, align = "right", fill = NA), running_sd3(x, 
        wins), running_skew4(x, wins), running_kurt5(x, 
        wins), running_tstat(x, wins), running_zscored(x, 
        wins), running_sharpe(x, wins), running_apx_median(x, 
        wins), running_centered(x, wins), running_scaled(x, 
        wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr   min    lq  mean median    uq    max neval       cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32237 32922 34848  34867 35800  44407   100        h 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64806 66702 71894  68078 69204 145095   100         i
##                                           running_sum(x, wins)    70    77   104     80    87   1964   100 a        
##                                          running_mean(x, wins)    71    77    85     80    94    131   100 a        
##                                       roll::roll_sum(xm, wins)  1886  1907  2010   1962  2028   3916   100  b d     
##                                      roll::roll_mean(xm, wins)  2066  2103  2168   2137  2207   2618   100   cd     
##                                        roll::roll_sd(xm, wins)  5846  5903  6094   5973  6056   9136   100     e    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   350   372   511    398   447   2744   100 abc      
##                             RollingWindow::RollingSum(x, wins)   105   121   194    135   180   3088   100 a        
##                            RollingWindow::RollingMean(x, wins)   140   150   195    163   204   1900   100 a        
##                             RollingWindow::RollingStd(x, wins)   227   239   333    256   312   2511   100 ab       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1738  1923  2007   1958  2029   4244   100  b d     
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1728  1903  2102   2008  2066   4135   100  b d     
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10276 10724 11367  10903 11569  14919   100      f   
##                                           running_sd3(x, wins)   610   620   738    642   673   2975   100 abc      
##                                         running_skew4(x, wins)  3610  3645  3772   3698  3756   6118   100    d     
##                                         running_kurt5(x, wins)  5830  5884  6186   5995  6315   8585   100     e    
##                                         running_tstat(x, wins)   655   664   686    675   695    820   100 abc      
##                                       running_zscored(x, wins)   659   669   743    682   725   2972   100 abc      
##                                        running_sharpe(x, wins)   658   664   712    674   716   2418   100 abc      
##                                    running_apx_median(x, wins) 13783 13921 14433  14052 14501  18123   100       g  
##                                      running_centered(x, wins)   571   581   606    592   617    883   100 abc      
##                                        running_scaled(x, wins)   657   664   739    680   712   3037   100 abc
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

