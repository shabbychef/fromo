

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
##                                "d207edfb078f"                                      "x86_64"                                     "unknown" 
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
##         expr   min    lq  mean median    uq   max neval    cld
##       sum(x)    80    80    83     81    86   109   100 a     
##      mean(x)   162   164   173    166   177   263   100 a     
##        sd(x)   459   468   504    489   508   727   100  b    
##  skewness(x)  8271  8362  8874   8466  9178 11816   100    d  
##  kurtosis(x)  8089  8181  8653   8330  8971 11798   100    d  
##       sd3(x)   850   860   899    874   934  1070   100   c   
##     skew4(x)  8351  8415  8755   8532  8793 10774   100    d  
##     kurt5(x) 14992 15149 15462  15264 15554 18098   100     e 
##     dumbk(x) 17028 17246 18218  17703 18796 23825   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16086 16345 18051  18532 19380 24590   100   c
##                                                               sd3(x, wts = w)   981   999  1107   1099  1145  1473   100 a  
##                                                                 slow_sd(x, w)  1405  1477  2277   1916  2882  4918   100  b
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
##  as.centsums(x1, 1)  187  189  222    203  232  448   100  b  
##  as.centsums(x1, 2)  113  116  136    127  137  273   100 a   
##  as.centsums(x1, 3)  850  854  937    877  956 1596   100   c 
##  as.centsums(x1, 4) 1482 1489 1691   1572 1760 2907   100    d
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
##  c(obj1, obj2)  14 15   18     15 16 174   100   b
##  obj3 %-% obj1  11 12   14     12 13  95   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.2  2.6    2.3 2.5  26   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.2    2.0 2.2  14   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.4    2.1 2.2  18   100   a
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
##  as.centcosums(x1, max_ord)  54 56   61     57 59 147   100   b
##             mobj3 %-% mobj1  17 18   21     20 21  38   100  a
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
##                                                           expr   min    lq  mean median    uq    max neval       cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32385 36298 41034  39689 44312 141198   100        h 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64487 71008 79715  76664 84178 188827   100         i
##                                           running_sum(x, wins)    74    83   121     93   112   2067   100 a        
##                                          running_mean(x, wins)    74    84   120     90   108   2251   100 a        
##                                       roll::roll_sum(xm, wins)  1895  2000  2263   2194  2481   4008   100 a c      
##                                      roll::roll_mean(xm, wins)  2074  2183  2393   2311  2495   3764   100  bc      
##                                        roll::roll_sd(xm, wins)  5737  6075  6642   6538  7090   8828   100     e    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   367   413   569    468   567   2913   100 ab       
##                             RollingWindow::RollingSum(x, wins)   115   153   312    188   232   3256   100 ab       
##                            RollingWindow::RollingMean(x, wins)   150   196   351    233   280   3942   100 ab       
##                             RollingWindow::RollingStd(x, wins)   232   274   453    329   383   2947   100 ab       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1744  2025  2228   2172  2391   2961   100 a c      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1753  2024  2340   2232  2501   5346   100  bc      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10241 11621 13957  12879 14357 101950   100      f   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   381   411   456    436   491    616   100 ab       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   447   482   543    522   593    800   100 ab       
##                                           running_sd3(x, wins)   620   667   750    717   811   1126   100 ab       
##                                          running_skew(x, wins)  3573  3770  4177   4096  4530   5544   100   cd     
##                                         running_skew4(x, wins)  3637  3888  4326   4191  4612   6804   100   cd     
##                                          running_kurt(x, wins)  5222  5728  6247   6184  6610   8250   100    de    
##                                         running_kurt5(x, wins)  5864  6194  6914   6827  7363  10098   100     e    
##                                         running_tstat(x, wins)   670   713   821    764   841   2665   100 ab       
##                                       running_zscored(x, wins)   697   742   841    805   920   1303   100 ab       
##                                        running_sharpe(x, wins)   687   734   843    784   872   3686   100 ab       
##                                    running_apx_median(x, wins) 14074 14863 16673  16493 18132  25050   100       g  
##                                      running_centered(x, wins)   577   628   696    673   744    956   100 ab       
##                                        running_scaled(x, wins)   694   752   837    801   903   1322   100 ab
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
    ungroup() %>% mutate(is_numeraire = grepl("^sum\\(x\\)$", 
    expr)) %>% arrange(!is_numeraire) %>% group_by(sernum) %>% 
    mutate(numv = first(meantime)) %>% ungroup() %>% 
    mutate(normalized = meantime/numv) %>% arrange(sernum) %>% 
    group_by(expr) %>% mutate(first_norm = first(normalized)) %>% 
    ungroup() %>% mutate(relchange = normalized/first_norm)

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
allt %>% arrange(sernum) %>% group_by(expr) %>% mutate(perfo = mean(relchange, 
    na.rm = TRUE), first_mean = first(meantime), last_mean = last(meantime)) %>% 
    ungroup() %>% distinct(expr, .keep_all = TRUE) %>% 
    select(expr, first_mean, last_mean, perfo) %>% 
    mutate(rel_mean = last_mean/first_mean) %>% arrange(desc(perfo)) %>% 
    head(n = 10) %>% kable()
```



|expr                                                        | first_mean| last_mean| perfo| rel_mean|
|:-----------------------------------------------------------|----------:|---------:|-----:|--------:|
|running_sum(x, wins)                                        |    4.3e+04|   9.2e+04|  1.21|     2.15|
|running_skew(x, wins)                                       |    3.3e+06|   4.2e+06|  1.03|     1.26|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) |    1.2e+07|   1.3e+07|  1.01|     1.10|
|sum(x)                                                      |    4.8e+04|   9.9e+04|  1.00|     2.04|
|RollingWindow::RollingStd(x, wins)                          |    3.7e+05|   4.3e+05|  0.99|     1.17|
|running_sd3(x, wins)                                        |    7.2e+05|   7.9e+05|  0.99|     1.10|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) |    4.3e+05|   4.4e+05|  0.96|     1.03|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   |    5.2e+05|   5.3e+05|  0.96|     1.02|
|running_kurt(x, wins)                                       |    5.9e+06|   5.8e+06|  0.92|     0.99|
|running_sum(x, wins, robust = FALSE)                        |    5.3e+04|   3.5e+04|  0.77|     0.65|


