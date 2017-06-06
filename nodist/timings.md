

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
##                                "6a0d2658b42d"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4200      
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
##       sum(x)    80    81    90     83    94   160   100 a      
##      mean(x)   162   168   187    179   197   309   100 a      
##        sd(x)   459   491   541    518   574   841   100  b     
##  skewness(x)  8320  8743  9368   9137  9706 12446   100     e  
##  kurtosis(x)  8135  8561  9080   8908  9518 11857   100    d   
##       sd3(x)   856   894   961    922  1033  1378   100   c    
##     skew4(x)  8569  8893  9307   9203  9565 11600   100    de  
##     kurt5(x) 15483 16122 16707  16647 17038 19407   100      f 
##     dumbk(x) 17398 18121 19045  18593 19849 22935   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16442 16892 17336  17262 17603 19846   100   c
##                                                               sd3(x, wts = w)   981  1001  1057   1025  1073  1338   100 a  
##                                                                 slow_sd(x, w)  1407  1445  2110   1776  2774  4367   100  b
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
##  as.centsums(x1, 1)  189  192  217    198  230  321   100  b  
##  as.centsums(x1, 2)  114  117  132    122  134  301   100 a   
##  as.centsums(x1, 3)  863  880  941    905  972 1335   100   c 
##  as.centsums(x1, 4) 1506 1535 1645   1585 1704 2161   100    d
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
##  c(obj1, obj2)  15 16   23     17 23 210   100   b
##  obj3 %-% obj1  12 13   17     14 15 113   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  3.2    2.4 2.6 32.6   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.2  2.9    2.3 2.5 25.7   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.1  2.3    2.2 2.4  7.5   100   a
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
##  as.centcosums(x1, max_ord)  54 56   66     58 69 183   100   b
##             mobj3 %-% mobj1  17 19   22     20 22  52   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 33964 34873 37135  36403 37448 120817   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 68582 70782 74508  72189 73541 154859   100           k
##                                           running_sum(x, wins)    73    77    88     83    97    157   100 a          
##                                          running_mean(x, wins)    73    78   112     86    97   2272   100 a          
##                                       roll::roll_sum(xm, wins)  1790  1871  2061   1939  2113   4774   100  bcd       
##                                      roll::roll_mean(xm, wins)  1953  2031  2251   2150  2327   4311   100   c e      
##                                        roll::roll_sd(xm, wins)  5425  5577  5972   5892  6325   7021   100       g    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   358   400   601    446   507   2919   100 a c        
##                             RollingWindow::RollingSum(x, wins)   111   138  1123    181   215  84625   100 a c        
##                            RollingWindow::RollingMean(x, wins)   139   165   344    210   241   2632   100 ab         
##                             RollingWindow::RollingStd(x, wins)   227   258   444    297   347   2762   100 a c        
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1752  1941  2016   2000  2072   2711   100  bc        
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1761  1985  2059   2046  2112   2907   100  bcd       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10420 10942 12057  11717 12909  15863   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   356   364   400    384   418    618   100 a c        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   414   425   495    440   504   2674   100 a c        
##                                           running_sd3(x, wins)   594   603   675    624   684   2704   100 a c        
##                                          running_skew(x, wins)  3626  3725  3879   3820  3976   4930   100    de      
##                                         running_skew4(x, wins)  3648  3733  3982   3848  4076   6803   100     ef     
##                                          running_kurt(x, wins)  5313  5492  5748   5669  5903   7058   100      fg    
##                                         running_kurt5(x, wins)  5929  6091  6374   6326  6542   9136   100       g    
##                                         running_tstat(x, wins)   628   637   684    656   693   1040   100 a c        
##                                       running_zscored(x, wins)   631   645   683    661   697    987   100 a c        
##                                        running_sharpe(x, wins)   630   642   698    660   689   2721   100 a c        
##                                    running_apx_median(x, wins) 14064 14574 14997  14924 15259  17368   100         i  
##                                      running_centered(x, wins)   554   567   636    586   623   2676   100 a c        
##                                        running_scaled(x, wins)   627   639   718    654   682   3549   100 a c
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

allt <- data.frame(fname = dir(".", "*.csv"), stringsAsFactors = FALSE) %>% 
    dplyr::filter(grepl("^timings_\\d.+\\d+.csv$", 
        fname)) %>% group_by(fname) %>% mutate(tims = list(readr::read_csv(fname))) %>% 
    ungroup() %>% tidyr::unnest() %>% mutate(sernum = gsub("^timings_(.+).csv$", 
    "\\1", fname)) %>% dplyr::select(-fname) %>% group_by(sernum, 
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
|running_sum(x, wins)                                                         |       0.51|      1.10|        2.18|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    155.95|        1.11|
|running_apx_median(x, wins)                                                  |     179.28|    194.57|        1.09|
|kurt5(x)                                                                     |     173.12|    183.06|        1.06|
|skew4(x)                                                                     |      97.66|    102.73|        1.05|
|mobj3 %-% mobj1                                                              |       0.26|      0.27|        1.04|
|running_sd3(x, wins)                                                         |       8.72|      9.11|        1.04|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      6.51|        1.03|
|running_kurt5(x, wins)                                                       |      76.25|     78.50|        1.03|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      5.35|        1.03|
|running_skew4(x, wins)                                                       |      48.07|     49.40|        1.03|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      3.23|        1.03|
|running_scaled(x, wins)                                                      |       9.03|      9.27|        1.03|
|roll::roll_sd(xm, wins)                                                      |      77.75|     79.68|        1.02|
|running_sharpe(x, wins)                                                      |       9.22|      9.38|        1.02|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.76|        1.02|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     27.18|        1.01|
|mean(x)                                                                      |       2.07|      2.09|        1.01|
|skewness(x)                                                                  |     100.16|    100.82|        1.01|
|running_tstat(x, wins)                                                       |       9.32|      9.38|        1.01|
|kurtosis(x)                                                                  |      98.46|     98.79|        1.00|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.49|        1.00|
|sd3(x)                                                                       |      10.52|     10.53|        1.00|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|running_kurt(x, wins)                                                        |      70.61|     70.61|        1.00|
|running_skew(x, wins)                                                        |      39.66|     39.66|        1.00|
|roll::roll_mean(xm, wins)                                                    |      28.39|     28.35|        1.00|
|dumbk(x)                                                                     |     209.85|    208.23|        0.99|
|sd(x)                                                                        |       5.82|      5.75|        0.99|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        0.99|
|running_zscored(x, wins)                                                     |       9.49|      9.31|        0.98|
|running_centered(x, wins)                                                    |       7.68|      7.50|        0.98|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     25.45|        0.98|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        0.97|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    191.57|        0.95|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.56|        0.95|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    887.26|        0.94|
|roll::roll_sum(xm, wins)                                                     |      27.28|     25.50|        0.93|
|sd3(x, wts = w)                                                              |      13.06|     11.73|        0.90|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    449.24|        0.90|
|as.centsums(x1, 4)                                                           |      20.20|     18.05|        0.89|
|obj3 %-% obj1                                                                |       0.19|      0.16|        0.87|
|as.centsums(x1, 3)                                                           |      11.79|     10.28|        0.87|
|c(obj1, obj2)                                                                |       0.25|      0.21|        0.87|
|as.centsums(x1, 1)                                                           |       2.83|      2.29|        0.81|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      2.99|        0.80|
|as.centsums(x1, 2)                                                           |       1.83|      1.43|        0.78|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.77|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_sum(x, wins, robust = TRUE)                                          |       1.88|      0.90|        0.48|


