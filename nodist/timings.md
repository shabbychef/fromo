

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
##                                "4a160b287f08"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4500      
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
##         expr   min    lq  mean median    uq   max neval    cld
##       sum(x)    80    82    91     83    94   219   100 a     
##      mean(x)   162   172   199    189   220   346   100 ab    
##        sd(x)   461   498   539    520   557   929   100  b    
##  skewness(x)  8377  8953  9591   9309  9958 14205   100    d  
##  kurtosis(x)  8205  8645  9294   9037  9601 14278   100    d  
##       sd3(x)   849   892   969    929  1028  1860   100   c   
##     skew4(x)  8663  9016  9575   9410  9833 13855   100    d  
##     kurt5(x) 15767 16568 17149  16839 17299 25018   100     e 
##     dumbk(x) 17458 18355 19454  19013 20283 26455   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16587 16912 17603  17426 18160 19969   100   c
##                                                               sd3(x, wts = w)   982  1008  1068   1041  1076  1385   100 a  
##                                                                 slow_sd(x, w)  1402  1488  2133   1721  2749  4518   100  b
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
##  as.centsums(x1, 1)  191  196  215    198  210  542   100  b  
##  as.centsums(x1, 2)  114  115  129    118  128  302   100 a   
##  as.centsums(x1, 3)  857  866  902    886  916 1050   100   c 
##  as.centsums(x1, 4) 1503 1528 1583   1551 1592 3037   100    d
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
##  c(obj1, obj2)  14 15   19     16 16 197   100   b
##  obj3 %-% obj1  11 12   15     13 13  96   100  a
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
##    join_cent_sums(rs1, rs2) 2.3 2.5  3.0    2.6 2.8  27   100   a
##  unjoin_cent_sums(rs3, rs2) 2.1 2.2  2.5    2.3 2.5  16   100   a
##  unjoin_cent_sums(rs3, rs1) 2.1 2.2  2.6    2.3 2.5  25   100   a
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
##  as.centcosums(x1, max_ord)  52 55   60     56 65 149   100   b
##             mobj3 %-% mobj1  17 18   21     19 21  46   100  a
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
##                                                           expr   min    lq  mean median    uq    max neval        cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32383 33074 35343  34908 35363 109327   100         i 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 65059 66156 70273  68181 68999 148598   100          j
##                                           running_sum(x, wins)   104   109   117    112   119    179   100 a         
##                                          running_mean(x, wins)   106   109   143    113   122   2591   100 a         
##                                       roll::roll_sum(xm, wins)  1789  1845  1916   1881  1943   2493   100  bc       
##                                      roll::roll_mean(xm, wins)  1966  1997  2092   2041  2109   4248   100   cd      
##                                        roll::roll_sd(xm, wins)  5450  5512  5847   5589  5810  16484   100      f    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   361   386   613    421   453   2722   100 a c       
##                             RollingWindow::RollingSum(x, wins)   109   135  1000    146   190  77222   100 a c       
##                            RollingWindow::RollingMean(x, wins)   142   163   260    182   230   2459   100 ab        
##                             RollingWindow::RollingStd(x, wins)   226   258   396    277   324   2613   100 a c       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1756  1832  1957   1965  2013   4198   100  bc       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1737  1821  1997   2007  2062   4075   100  bc       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10409 10794 11625  11048 12712  14635   100       g   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   355   361   402    371   390   2271   100 a c       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   414   420   439    426   444    731   100 a c       
##                                           running_sd3(x, wins)   586   592   618    606   619    971   100 a c       
##                                          running_skew(x, wins)  3642  3682  3793   3734  3806   6331   100    de     
##                                         running_skew4(x, wins)  3639  3674  3805   3727  3811   6060   100    de     
##                                          running_kurt(x, wins)  5339  5368  5499   5431  5512   6681   100     ef    
##                                         running_kurt5(x, wins)  5935  6025  6224   6093  6210   8791   100      f    
##                                         running_tstat(x, wins)   623   630   679    643   669   2368   100 a c       
##                                       running_zscored(x, wins)   622   634   665    652   673   1057   100 a c       
##                                        running_sharpe(x, wins)   620   626   662    639   668   1068   100 a c       
##                                    running_apx_median(x, wins) 13979 14120 14529  14268 14554  17211   100        h  
##                                      running_centered(x, wins)   548   558   584    574   592    770   100 a c       
##                                        running_scaled(x, wins)   621   627   648    635   648    829   100 a c
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
|RollingWindow::RollingSum(x, wins)                                           |       3.72|     10.97|        2.95|
|running_sum(x, wins)                                                         |       0.51|      1.28|        2.53|
|kurt5(x)                                                                     |     173.12|    188.17|        1.09|
|skew4(x)                                                                     |      97.66|    105.06|        1.08|
|mean(x)                                                                      |       2.07|      2.19|        1.06|
|skewness(x)                                                                  |     100.16|    105.23|        1.05|
|running_skew(x, wins)                                                        |      39.66|     41.62|        1.05|
|kurtosis(x)                                                                  |      98.46|    101.98|        1.04|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        1.03|
|dumbk(x)                                                                     |     209.85|    213.46|        1.02|
|sd(x)                                                                        |       5.82|      5.92|        1.02|
|sd3(x)                                                                       |      10.52|     10.63|        1.01|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        0.99|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.35|        0.97|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.72|        0.97|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    193.15|        0.96|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    127.56|        0.91|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      2.85|        0.91|
|sd3(x, wts = w)                                                              |      13.06|     11.72|        0.90|
|running_kurt5(x, wins)                                                       |      76.25|     68.30|        0.90|
|running_apx_median(x, wins)                                                  |     179.28|    159.42|        0.89|
|mobj3 %-% mobj1                                                              |       0.26|      0.23|        0.88|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.66|        0.87|
|running_skew4(x, wins)                                                       |      48.07|     41.75|        0.87|
|as.centsums(x1, 4)                                                           |      20.20|     17.37|        0.86|
|running_kurt(x, wins)                                                        |      70.61|     60.34|        0.85|
|c(obj1, obj2)                                                                |       0.25|      0.21|        0.85|
|obj3 %-% obj1                                                                |       0.19|      0.16|        0.85|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      4.41|        0.85|
|as.centsums(x1, 3)                                                           |      11.79|      9.90|        0.84|
|running_centered(x, wins)                                                    |       7.68|      6.41|        0.84|
|as.centsums(x1, 1)                                                           |       2.83|      2.36|        0.83|
|roll::roll_sd(xm, wins)                                                      |      77.75|     64.15|        0.83|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     21.48|        0.82|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.82|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     21.91|        0.82|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    771.09|        0.82|
|roll::roll_mean(xm, wins)                                                    |      28.39|     22.96|        0.81|
|running_tstat(x, wins)                                                       |       9.32|      7.45|        0.80|
|running_sharpe(x, wins)                                                      |       9.22|      7.26|        0.79|
|running_scaled(x, wins)                                                      |       9.03|      7.11|        0.79|
|running_sd3(x, wins)                                                         |       8.72|      6.78|        0.78|
|as.centsums(x1, 2)                                                           |       1.83|      1.42|        0.77|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    387.81|        0.77|
|roll::roll_sum(xm, wins)                                                     |      27.28|     21.02|        0.77|
|running_zscored(x, wins)                                                     |       9.49|      7.29|        0.77|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      4.81|        0.76|
|running_mean(x, wins)                                                        |       2.36|      1.57|        0.66|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|


