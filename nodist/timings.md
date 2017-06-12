

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
##                                "28f6ec8a4353"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4510      
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
##       sum(x)    80    80    87     81    89   155   100 a     
##      mean(x)   162   164   175    167   175   338   100 ab    
##        sd(x)   460   466   507    486   525   855   100  b    
##  skewness(x)  8281  8376  9011   8630  9201 12163   100    d  
##  kurtosis(x)  8109  8209  8828   8480  9173 11952   100    d  
##       sd3(x)   849   856   909    876   941  1083   100   c   
##     skew4(x)  8540  8616  9128   8813  9432 12141   100    d  
##     kurt5(x) 15406 15565 16286  15832 16509 21590   100     e 
##     dumbk(x) 17062 17446 18644  18067 19574 23281   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16268 16590 18613  18033 19720 28964   100   c
##                                                               sd3(x, wts = w)   981  1000  1109   1082  1149  1568   100 a  
##                                                                 slow_sd(x, w)  1415  1563  2374   2030  3236  5290   100  b
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
##  as.centsums(x1, 1)  191  197  209    198  223  296   100  b  
##  as.centsums(x1, 2)  114  115  125    117  130  320   100 a   
##  as.centsums(x1, 3)  859  864  914    870  924 1216   100   c 
##  as.centsums(x1, 4) 1504 1507 1585   1523 1559 2153   100    d
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
##  c(obj1, obj2)  14 15   18     16 16 178   100   b
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
##    join_cent_sums(rs1, rs2) 2.3 2.4  2.9    2.5 2.7  27   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.2  2.5    2.2 2.5  17   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.2  2.5    2.3 2.5  14   100   a
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
##  as.centcosums(x1, max_ord)  55 59   75     66 78 228   100   b
##             mobj3 %-% mobj1  18 21   26     23 27  70   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32123 33664 38030  35450 40129 111857   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64396 66983 74071  69176 76926 159247   100           k
##                                           running_sum(x, wins)   104   108   121    115   129    202   100 a          
##                                          running_mean(x, wins)   105   109   148    118   130   2303   100 ab         
##                                       roll::roll_sum(xm, wins)  1862  1938  2138   2009  2172   4060   100  bc e      
##                                      roll::roll_mean(xm, wins)  2026  2124  2384   2209  2450   6452   100    de      
##                                        roll::roll_sd(xm, wins)  5676  5796  6174   5883  6527   8151   100       g    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   359   386   654    433   568   3053   100 a cd       
##                             RollingWindow::RollingSum(x, wins)   110   131  1019    160   205  76877   100 a cd       
##                            RollingWindow::RollingMean(x, wins)   141   157   279    197   245   2571   100 a c        
##                             RollingWindow::RollingStd(x, wins)   225   251   395    286   328   2586   100 a cd       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1750  1945  2150   2020  2257   4974   100   c e      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1730  1977  2120   2038  2255   4136   100  bc e      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10368 11042 12463  12605 13183  20783   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   355   361   417    390   413   2248   100 a cd       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   410   425   481    452   489    855   100 a cd       
##                                           running_sd3(x, wins)   586   597   668    638   682   1276   100 a cd       
##                                          running_skew(x, wins)  3631  3684  4091   3831  4257   5866   100     ef     
##                                         running_skew4(x, wins)  3639  3683  4047   3789  4091   7580   100     ef     
##                                          running_kurt(x, wins)  5320  5367  5889   5473  6215   9933   100      fg    
##                                         running_kurt5(x, wins)  5949  6017  6507   6154  6626  10240   100       g    
##                                         running_tstat(x, wins)   620   638   731    676   720   3376   100 a cd       
##                                       running_zscored(x, wins)   625   631   698    659   709   1250   100 a cd       
##                                        running_sharpe(x, wins)   620   628   706    674   716   1274   100 a cd       
##                                    running_apx_median(x, wins) 13901 14181 15592  15034 16280  24211   100         i  
##                                      running_centered(x, wins)   549   555   608    589   626    885   100 a cd       
##                                        running_scaled(x, wins)   621   629   694    657   713   1201   100 a cd
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
|RollingWindow::RollingSum(x, wins)                                           |       3.72|     11.66|        3.14|
|running_sum(x, wins)                                                         |       0.51|      1.39|        2.74|
|running_skew(x, wins)                                                        |      39.66|     46.83|        1.18|
|mobj3 %-% mobj1                                                              |       0.26|      0.30|        1.17|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.85|        1.13|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      7.48|        1.08|
|kurt5(x)                                                                     |     173.12|    186.41|        1.08|
|skew4(x)                                                                     |      97.66|    104.48|        1.07|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    213.04|        1.06|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        1.05|
|skewness(x)                                                                  |     100.16|    103.15|        1.03|
|kurtosis(x)                                                                  |      98.46|    101.05|        1.03|
|dumbk(x)                                                                     |     209.85|    213.40|        1.02|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      3.19|        1.02|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    142.65|        1.01|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.52|        1.01|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|sd(x)                                                                        |       5.82|      5.81|        1.00|
|running_apx_median(x, wins)                                                  |     179.28|    178.47|        1.00|
|sd3(x)                                                                       |      10.52|     10.41|        0.99|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        0.98|
|running_kurt5(x, wins)                                                       |      76.25|     74.48|        0.98|
|sd3(x, wts = w)                                                              |      13.06|     12.69|        0.97|
|mean(x)                                                                      |       2.07|      2.00|        0.97|
|running_skew4(x, wins)                                                       |      48.07|     46.32|        0.96|
|roll::roll_mean(xm, wins)                                                    |      28.39|     27.29|        0.96|
|running_kurt(x, wins)                                                        |      70.61|     67.40|        0.95|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     24.61|        0.94|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      4.78|        0.92|
|roll::roll_sd(xm, wins)                                                      |      77.75|     70.67|        0.91|
|running_centered(x, wins)                                                    |       7.68|      6.95|        0.91|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     24.26|        0.91|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    847.82|        0.90|
|as.centsums(x1, 4)                                                           |      20.20|     18.14|        0.90|
|running_tstat(x, wins)                                                       |       9.32|      8.36|        0.90|
|roll::roll_sum(xm, wins)                                                     |      27.28|     24.47|        0.90|
|as.centsums(x1, 3)                                                           |      11.79|     10.46|        0.89|
|running_scaled(x, wins)                                                      |       9.03|      7.94|        0.88|
|running_sd3(x, wins)                                                         |       8.72|      7.65|        0.88|
|running_sharpe(x, wins)                                                      |       9.22|      8.08|        0.88|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      5.51|        0.87|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    435.29|        0.87|
|obj3 %-% obj1                                                                |       0.19|      0.16|        0.85|
|c(obj1, obj2)                                                                |       0.25|      0.21|        0.85|
|as.centsums(x1, 1)                                                           |       2.83|      2.40|        0.85|
|running_zscored(x, wins)                                                     |       9.49|      7.99|        0.84|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.84|
|as.centsums(x1, 2)                                                           |       1.83|      1.43|        0.78|
|running_mean(x, wins)                                                        |       2.36|      1.69|        0.72|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|


