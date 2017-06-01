

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
##                                "9e78e42a50b8"                                      "x86_64"                                     "unknown" 
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
##       sum(x)    80    80    86     82    87   160   100 a      
##      mean(x)   162   165   178    168   182   295   100 a      
##        sd(x)   459   473   498    485   502   776   100  b     
##  skewness(x)  8324  8515  9031   8729  9157 11830   100     e  
##  kurtosis(x)  8141  8343  8684   8484  8743 11209   100    d   
##       sd3(x)   855   872   933    895   995  1209   100   c    
##     skew4(x)  8375  8517  8784   8671  8890 10186   100    de  
##     kurt5(x) 15111 15374 15831  15553 15961 20315   100      f 
##     dumbk(x) 17328 17710 18607  18195 19241 23145   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16255 16640 17002  16973 17214 20016   100   c
##                                                               sd3(x, wts = w)   981   998  1044   1014  1067  1316   100 a  
##                                                                 slow_sd(x, w)  1402  1443  2015   1557  2696  4337   100  b
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
##  as.centsums(x1, 1)  187  188  196    190  201  295   100  b  
##  as.centsums(x1, 2)  113  115  123    116  125  277   100 a   
##  as.centsums(x1, 3)  845  852  879    870  879 1084   100   c 
##  as.centsums(x1, 4) 1473 1496 1532   1514 1546 2024   100    d
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
##  c(obj1, obj2)  14 15   19     16 16 192   100   b
##  obj3 %-% obj1  11 12   14     12 13 102   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.2  2.7    2.3 2.5  26   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.3    2.1 2.3  15   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.3    2.1 2.2  14   100   a
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
##  as.centcosums(x1, max_ord)  51 52   58     53 63 140   100   b
##             mobj3 %-% mobj1  16 18   20     19 20  32   100  a
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
##                                                           expr   min    lq  mean median    uq    max neval      cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32271 32782 35441  34526 35290 105955   100       g 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 63652 65343 71540  67036 70076 193077   100        h
##                                           running_sum(x, wins)    73    77   109     84    97   2121   100 a       
##                                          running_mean(x, wins)    73    77    85     79    92    144   100 a       
##                                       roll::roll_sum(xm, wins)  1788  1828  1998   1917  1996   4781   100 a c     
##                                      roll::roll_mean(xm, wins)  1958  2031  2244   2107  2242   6014   100  bc     
##                                        roll::roll_sd(xm, wins)  5419  5532  5903   5646  5963   7735   100    de   
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   353   373   466    405   464   2207   100 ab      
##                             RollingWindow::RollingSum(x, wins)   108   120   254    146   203   2488   100 ab      
##                            RollingWindow::RollingMean(x, wins)   138   156   237    181   231   2289   100 ab      
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1739  1927  2044   1985  2071   4182   100 a c     
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1724  1926  2102   2026  2098   4016   100 a c     
##                                         running_skew4(x, wins)  3613  3708  3963   3766  3896   9024   100   cd    
##                                         running_kurt5(x, wins)  5836  5978  6337   6102  6424   9308   100     e   
##                                         running_tstat(x, wins)   658   665   761    682   773   2873   100 ab      
##                                       running_zscored(x, wins)   662   670   725    693   756   1140   100 ab      
##                                        running_sharpe(x, wins)   657   672   720    684   761   1005   100 ab      
##                                    running_apx_median(x, wins) 14034 14202 14929  14396 15084  20009   100      f  
##                                      running_centered(x, wins)   572   580   631    601   645   1043   100 ab      
##                                        running_scaled(x, wins)   658   669   734    685   738   2540   100 ab
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

