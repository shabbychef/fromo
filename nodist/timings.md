

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
##                                "b70aaf3783fc"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.5000      
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
##       sum(x)    80    87    93     90    98   141   100 a      
##      mean(x)   168   183   198    195   208   314   100 a      
##        sd(x)   489   524   554    538   570   778   100  b     
##  skewness(x)  8522  9123  9610   9359  9913 13839   100    d   
##  kurtosis(x)  8191  8963  9474   9162  9716 14485   100    d   
##       sd3(x)   862   941   981    955  1015  1224   100   c    
##     skew4(x)  8850  9651  9930   9780 10056 11819   100     e  
##     kurt5(x) 16020 16976 17486  17222 17585 21824   100      f 
##     dumbk(x) 17403 18958 19935  19343 20803 25003   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16495 16622 16981  16725 17062 19797   100   c
##                                                               sd3(x, wts = w)   980   983  1010    991  1014  1172   100 a  
##                                                                 slow_sd(x, w)  1389  1424  2009   1494  2788  3977   100  b
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
##  as.centsums(x1, 1)  193  194  203    195  209  258   100  b  
##  as.centsums(x1, 2)  114  115  123    117  126  309   100 a   
##  as.centsums(x1, 3)  881  886  915    899  944 1003   100   c 
##  as.centsums(x1, 4) 1529 1537 1583   1557 1599 1858   100    d
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
##  c(obj1, obj2)  14 15   19     15 16 193   100   b
##  obj3 %-% obj1  11 12   14     12 13  93   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  2.8    2.4 2.6  25   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.3    2.1 2.3  15   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.3    2.1 2.3  14   100   a
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
##  as.centcosums(x1, max_ord)  51 55   60     56 64 154   100   b
##             mobj3 %-% mobj1  16 18   20     19 20  41   100  a
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
        n = wins, align = "right", fill = NA), running_sd(x, 
        wins, na_rm = FALSE, restart_period = 50000L), 
    running_sd(x, wins, na_rm = TRUE, restart_period = 1000L), 
    running_sd3(x, wins), running_skew(x, wins), running_skew4(x, 
        wins), running_kurt(x, wins), running_kurt5(x, 
        wins), running_tstat(x, wins), running_zscored(x, 
        wins), running_sharpe(x, wins), running_apx_median(x, 
        wins), running_centered(x, wins), running_scaled(x, 
        wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr   min    lq  mean median    uq    max neval         cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32568 33038 35422  34956 35765 109978   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64438 65544 70321  67541 68981 149314   100           k
##                                           running_sum(x, wins)   105   108   115    111   120    165   100 a          
##                                          running_mean(x, wins)   104   108   115    111   119    153   100 a          
##                                       roll::roll_sum(xm, wins)  1897  1937  2046   1977  2037   4224   100   cd       
##                                      roll::roll_mean(xm, wins)  2064  2103  2225   2161  2219   4403   100    de      
##                                        roll::roll_sd(xm, wins)  5719  5785  5981   5849  5906   9198   100       g    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   367   391   609    410   450   2984   100 a  d       
##                             RollingWindow::RollingSum(x, wins)   118   137   995    156   185  78855   100 a  d       
##                            RollingWindow::RollingMean(x, wins)   146   167   259    183   225   2489   100 ab         
##                             RollingWindow::RollingStd(x, wins)   230   255   417    270   319   2743   100 abc        
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1747  1840  1961   1957  2011   2379   100  b d       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1734  1982  1999   2025  2053   2359   100  b d       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10298 10691 11579  11170 12652  14668   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   356   362   401    370   392   2525   100 abc        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   416   423   441    430   452    651   100 abc        
##                                           running_sd3(x, wins)   597   607   652    622   647   2815   100 a  d       
##                                          running_skew(x, wins)  3680  3698  3837   3752  3845   5757   100     e      
##                                         running_skew4(x, wins)  3693  3731  3869   3768  3909   6130   100     ef     
##                                          running_kurt(x, wins)  5389  5443  5649   5502  5743   7684   100      fg    
##                                         running_kurt5(x, wins)  6012  6053  6299   6117  6393   8775   100       g    
##                                         running_tstat(x, wins)   628   637   662    648   676    788   100 a  d       
##                                       running_zscored(x, wins)   631   639   689    651   676   2990   100 a  d       
##                                        running_sharpe(x, wins)   626   633   663    647   682    829   100 a  d       
##                                    running_apx_median(x, wins) 14042 14167 14641  14261 14838  17958   100         i  
##                                      running_centered(x, wins)   587   594   620    604   637    748   100 a  d       
##                                        running_scaled(x, wins)   627   636   669    654   686    850   100 a  d
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

## performance regressions?

Can we see them here? load all the timing data to check.


```r
library(magrittr)
library(dplyr)
library(readr)
library(tidyr)
library(knitr)

mysernum <- as.character(packageVersion("fromo"))
allt <- data.frame(fname = dir(".", "*.csv"), stringsAsFactors = FALSE) %>% 
    dplyr::filter(grepl("^timings_\\d.+\\d+.csv$", 
        fname)) %>% group_by(fname) %>% mutate(tims = list(readr::read_csv(fname))) %>% 
    ungroup() %>% tidyr::unnest() %>% mutate(sernum = gsub("^timings_(.+).csv$", 
    "\\1", fname)) %>% dplyr::select(-fname) %>% rbind(resdf %>% 
    dplyr::mutate(sernum = mysernum)) %>% group_by(sernum, 
    expr) %>% summarize(meantime = mean(time, na.rm = TRUE)) %>% 
    ungroup() %>% group_by(sernum) %>% mutate(sumx_time = median(ifelse(grepl("^sum\\(x\\)$", 
    expr), meantime, NA), na.rm = TRUE)) %>% ungroup() %>% 
    mutate(normalized = meantime/sumx_time) %>% arrange(sernum) %>% 
    group_by(expr) %>% mutate(first_norm = first(normalized), 
    last_norm = last(normalized)) %>% ungroup() %>% 
    mutate(relchange = normalized/first_norm, last_status = last_norm/first_norm)

library(ggplot2)
ph <- allt %>% ggplot(aes(sernum, normalized, group = expr, 
    color = expr)) + geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + labs(x = "release", y = "mean time taken, relative to sum(x)", 
    title = "fromo microbenchmark timings")
print(ph)
```

<img src="figure/all_timing_stats-1.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="600px" height="500px" />

```r
ph <- allt %>% ggplot(aes(sernum, relchange, group = expr, 
    color = expr)) + geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + labs(x = "release", y = "normalized time taken, relative to first iteration", 
    title = "fromo microbenchmark timings")
print(ph)
```

<img src="figure/all_timing_stats-2.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="600px" height="500px" />

```r
allt %>% select(-sernum, -relchange, -meantime, -sumx_time, 
    -normalized) %>% distinct(expr, .keep_all = TRUE) %>% 
    arrange(desc(last_status)) %>% head(n = 50) %>% 
    kable()
```



|expr                                                                         | first_norm| last_norm| last_status|
|:----------------------------------------------------------------------------|----------:|---------:|-----------:|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|     10.68|        2.87|
|running_sum(x, wins)                                                         |       0.51|      1.24|        2.45|
|skew4(x)                                                                     |      97.66|    106.66|        1.09|
|kurt5(x)                                                                     |     173.12|    187.82|        1.08|
|running_skew(x, wins)                                                        |      39.66|     41.21|        1.04|
|kurtosis(x)                                                                  |      98.46|    101.76|        1.03|
|skewness(x)                                                                  |     100.16|    103.22|        1.03|
|mean(x)                                                                      |       2.07|      2.13|        1.03|
|sd(x)                                                                        |       5.82|      5.95|        1.02|
|dumbk(x)                                                                     |     209.85|    214.12|        1.02|
|sd3(x)                                                                       |      10.52|     10.54|        1.00|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.48|        1.00|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.54|        0.94|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        0.93|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    182.39|        0.91|
|running_kurt5(x, wins)                                                       |      76.25|     67.66|        0.89|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    124.37|        0.88|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      2.78|        0.88|
|running_apx_median(x, wins)                                                  |     179.28|    157.26|        0.88|
|running_centered(x, wins)                                                    |       7.68|      6.66|        0.87|
|running_skew4(x, wins)                                                       |      48.07|     41.55|        0.86|
|running_kurt(x, wins)                                                        |      70.61|     60.67|        0.86|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.64|        0.85|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.02|        0.85|
|as.centsums(x1, 4)                                                           |      20.20|     17.01|        0.84|
|roll::roll_mean(xm, wins)                                                    |      28.39|     23.90|        0.84|
|as.centsums(x1, 3)                                                           |      11.79|      9.83|        0.83|
|sd3(x, wts = w)                                                              |      13.06|     10.85|        0.83|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      4.31|        0.83|
|mobj3 %-% mobj1                                                              |       0.26|      0.21|        0.83|
|roll::roll_sd(xm, wins)                                                      |      77.75|     64.24|        0.83|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     21.06|        0.81|
|c(obj1, obj2)                                                                |       0.25|      0.20|        0.81|
|roll::roll_sum(xm, wins)                                                     |      27.28|     21.98|        0.81|
|running_sd3(x, wins)                                                         |       8.72|      7.00|        0.80|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     21.47|        0.80|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    755.31|        0.80|
|running_scaled(x, wins)                                                      |       9.03|      7.18|        0.79|
|obj3 %-% obj1                                                                |       0.19|      0.15|        0.79|
|running_zscored(x, wins)                                                     |       9.49|      7.40|        0.78|
|running_sharpe(x, wins)                                                      |       9.22|      7.12|        0.77|
|as.centsums(x1, 1)                                                           |       2.83|      2.18|        0.77|
|running_tstat(x, wins)                                                       |       9.32|      7.11|        0.76|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    380.47|        0.76|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      4.74|        0.75|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.02|        0.73|
|as.centsums(x1, 2)                                                           |       1.83|      1.32|        0.72|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.24|        0.52|


