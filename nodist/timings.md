

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
##                                "9d552a41626c"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.6000      
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
##       sum(x)    80    80    87     81    88   161   100 a      
##      mean(x)   162   164   176    166   179   320   100 ab     
##        sd(x)   459   465   502    484   520   769   100  b     
##  skewness(x)  8293  8345  8906   8430  9006 12650   100    d   
##  kurtosis(x)  8105  8175  8636   8285  8794 11015   100    d   
##       sd3(x)  1463  1472  1561   1505  1599  2867   100   c    
##     skew4(x)  8831  8902  9613   9043  9952 16376   100     e  
##     kurt5(x) 15789 15951 16627  16131 16737 21468   100      f 
##     dumbk(x) 17062 17228 18366  17609 19036 22752   100       g
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
## Unit: milliseconds
##                                                                          expr  min   lq mean median   uq  max neval cld
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16.3 16.4 16.8   16.5 16.8 20.1   100   c
##                                                               sd3(x, wts = w)  1.6  1.6  1.7    1.7  1.7  2.2   100 a  
##                                                                 slow_sd(x, w)  1.4  1.4  2.0    1.5  2.7  4.4   100  b
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
##                expr  min   lq mean median   uq  max neval cld
##  as.centsums(x1, 1)  194  197  212    208  221  265   100 a  
##  as.centsums(x1, 2)  197  200  216    212  217  458   100 a  
##  as.centsums(x1, 3)  898  904  970    976  995 1155   100  b 
##  as.centsums(x1, 4) 1550 1556 1660   1656 1704 1978   100   c
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
##  c(obj1, obj2)  14 15   20     16 17 190   100   b
##  obj3 %-% obj1  11 12   14     13 13  96   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.4  2.8    2.5 2.7  26   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.2  3.0    2.3 2.4  39   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.2  2.6    2.3 2.5  25   100   a
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
##  as.centcosums(x1, max_ord)  50 52   57     54 61 150   100   b
##             mobj3 %-% mobj1  17 18   20     19 20  57   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32348 33173 36702  34692 35335 127263   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 63717 65329 70182  67089 70223 154537   100           k
##                                           running_sum(x, wins)   104   109   157    112   125   2128   100 ab         
##                                          running_mean(x, wins)   103   108   119    114   124    186   100 a          
##                                       roll::roll_sum(xm, wins)  1884  1913  2008   1943  1998   3970   100 abc        
##                                      roll::roll_mean(xm, wins)  2072  2101  2207   2127  2214   3232   100   cd       
##                                        roll::roll_sd(xm, wins)  5707  5760  6009   5829  5940   8137   100      fg    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   356   378   545    409   451   2706   100 abc        
##                             RollingWindow::RollingSum(x, wins)   111   127   215    145   186   2157   100 ab         
##                            RollingWindow::RollingMean(x, wins)   141   157   233    178   222   2299   100 ab         
##                             RollingWindow::RollingStd(x, wins)   225   245   305    263   311   2541   100 abc        
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1743  1799  1974   1947  2003   4037   100 abc        
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1721  1818  2055   2011  2084   4209   100  b d       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10126 10618 12273  10943 12710  85925   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   427   434   463    445   473    674   100 abc        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   560   586   617    598   622    972   100 abc        
##                                           running_sd3(x, wins)  1341  1384  1508   1425  1497   3493   100 abc        
##                                          running_skew(x, wins)  3704  3734  3933   3772  3896   7481   100    de      
##                                         running_skew4(x, wins)  4434  4501  4810   4591  4811   7425   100     ef     
##                                          running_kurt(x, wins)  5444  5475  5818   5550  5937   8380   100     e g    
##                                         running_kurt5(x, wins)  6795  6871  7194   6970  7100   9464   100       g    
##                                         running_tstat(x, wins)   893   899   977    925   968   2548   100 abc        
##                                       running_zscored(x, wins)   861   867   948    881   940   3151   100 abc        
##                                        running_sharpe(x, wins)   856   860   916    869   905   2903   100 abc        
##                                    running_apx_median(x, wins) 14574 14714 15553  14855 16121  22819   100         i  
##                                      running_centered(x, wins)   619   626   679    635   664   2601   100 abc        
##                                        running_scaled(x, wins)   858   862   920    878   928   1529   100 abc
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
|running_sum(x, wins)                                                         |       0.51|      1.81|        3.59|
|running_sd3(x, wins)                                                         |       8.72|     17.37|        1.99|
|sd3(x)                                                                       |      10.52|     17.98|        1.71|
|sd3(x, wts = w)                                                              |      13.06|     19.77|        1.51|
|as.centsums(x1, 2)                                                           |       1.83|      2.49|        1.36|
|running_tstat(x, wins)                                                       |       9.32|     11.25|        1.21|
|running_scaled(x, wins)                                                      |       9.03|     10.59|        1.17|
|running_skew4(x, wins)                                                       |      48.07|     55.41|        1.15|
|running_zscored(x, wins)                                                     |       9.49|     10.92|        1.15|
|running_sharpe(x, wins)                                                      |       9.22|     10.55|        1.14|
|running_skew(x, wins)                                                        |      39.66|     45.31|        1.14|
|skew4(x)                                                                     |      97.66|    110.75|        1.13|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      7.11|        1.13|
|kurt5(x)                                                                     |     173.12|    191.56|        1.11|
|running_kurt5(x, wins)                                                       |      76.25|     82.88|        1.09|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        1.03|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        1.03|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      5.34|        1.03|
|skewness(x)                                                                  |     100.16|    102.60|        1.02|
|running_centered(x, wins)                                                    |       7.68|      7.83|        1.02|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        1.01|
|kurtosis(x)                                                                  |      98.46|     99.50|        1.01|
|dumbk(x)                                                                     |     209.85|    211.59|        1.01|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    141.39|        1.01|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|running_apx_median(x, wins)                                                  |     179.28|    179.19|        1.00|
|sd(x)                                                                        |       5.82|      5.79|        0.99|
|mean(x)                                                                      |       2.07|      2.02|        0.98|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    193.53|        0.96|
|running_kurt(x, wins)                                                        |      70.61|     67.02|        0.95|
|as.centsums(x1, 3)                                                           |      11.79|     11.17|        0.95|
|as.centsums(x1, 4)                                                           |      20.20|     19.13|        0.95|
|c(obj1, obj2)                                                                |       0.25|      0.23|        0.92|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.28|        0.91|
|mobj3 %-% mobj1                                                              |       0.26|      0.23|        0.90|
|roll::roll_mean(xm, wins)                                                    |      28.39|     25.42|        0.90|
|roll::roll_sd(xm, wins)                                                      |      77.75|     69.23|        0.89|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     23.67|        0.88|
|obj3 %-% obj1                                                                |       0.19|      0.17|        0.88|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.66|        0.88|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     22.74|        0.87|
|as.centsums(x1, 1)                                                           |       2.83|      2.44|        0.86|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    808.56|        0.86|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      2.68|        0.85|
|roll::roll_sum(xm, wins)                                                     |      27.28|     23.13|        0.85|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    422.84|        0.84|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      3.51|        0.78|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      2.48|        0.67|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.37|        0.58|


