

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
##                                "4e4828ca9c20"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.5020      
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
##       sum(x)    79    80    82     80    81   124   100 a      
##      mean(x)   162   163   170    164   173   304   100 ab     
##        sd(x)   459   463   491    472   498   728   100  b     
##  skewness(x)  8280  8351  8951   8475  9070 14532   100    de  
##  kurtosis(x)  8103  8173  8655   8278  8821 11785   100    d   
##       sd3(x)  1474  1479  1558   1500  1573  2907   100   c    
##     skew4(x)  8836  8893  9186   8944  9118 12614   100     e  
##     kurt5(x) 15806 15904 16559  16043 16929 23117   100      f 
##     dumbk(x) 17008 17145 17944  17396 18701 23706   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16.2 16.3 16.7   16.4 16.6 20.0   100   c
##                                                               sd3(x, wts = w)  1.6  1.6  1.7    1.7  1.7  1.8   100 a  
##                                                                 slow_sd(x, w)  1.4  1.4  1.9    1.4  2.7  3.7   100  b
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
##  as.centsums(x1, 1)  188  190  197    191  205  238   100 a  
##  as.centsums(x1, 2)  199  200  209    201  216  396   100 a  
##  as.centsums(x1, 3)  895  898  924    900  972 1001   100  b 
##  as.centsums(x1, 4) 1545 1549 1586   1553 1576 1818   100   c
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
##  c(obj1, obj2)  14 15   18     16 16 180   100   a
##  obj3 %-% obj1  11 12   15     13 13  93   100   a
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
##    join_cent_sums(rs1, rs2) 2.5 2.6  3.0    2.7 2.9  27   100   a
##  unjoin_cent_sums(rs3, rs2) 2.3 2.5  2.7    2.6 2.7  16   100   a
##  unjoin_cent_sums(rs3, rs1) 2.3 2.4  2.9    2.5 2.7  21   100   a
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
##  as.centcosums(x1, max_ord)  64 66   72     67 71 163   100   b
##             mobj3 %-% mobj1  20 22   24     23 24  41   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32278 32952 36204  34718 35346 113530   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 63969 65735 68765  66866 68315 140801   100           k
##                                           running_sum(x, wins)   104   109   155    112   126   2084   100 a          
##                                          running_mean(x, wins)   103   108   116    111   118    169   100 a          
##                                       roll::roll_sum(xm, wins)  1886  1919  2023   1954  2027   4421   100  bc        
##                                      roll::roll_mean(xm, wins)  2066  2096  2163   2127  2186   2782   100   cd       
##                                        roll::roll_sd(xm, wins)  5705  5766  6021   5822  5910  13624   100      fg    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   364   378   533    403   443   2676   100 a c        
##                             RollingWindow::RollingSum(x, wins)   110   126   209    142   186   2163   100 a          
##                            RollingWindow::RollingMean(x, wins)   139   155   229    178   220   2302   100 a          
##                             RollingWindow::RollingStd(x, wins)   226   243   298    262   307   2511   100 ab         
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1740  1924  1967   1952  1996   4183   100  bc        
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1726  1891  2012   1997  2052   3842   100  bc        
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10129 10460 12007  10756 12388  84982   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   427   432   453    437   463    548   100 a c        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   522   541   563    552   578    676   100 a c        
##                                           running_sd3(x, wins)  1326  1356  1463   1383  1468   3489   100 a c        
##                                          running_skew(x, wins)  3708  3719  3857   3758  3802   6140   100    de      
##                                         running_skew4(x, wins)  4412  4494  4657   4563  4658   6498   100     ef     
##                                          running_kurt(x, wins)  5444  5466  5569   5489  5582   6521   100     e g    
##                                         running_kurt5(x, wins)  6760  6862  7116   6921  7058   9272   100       g    
##                                         running_tstat(x, wins)   891   898   947    911   935   2579   100 a c        
##                                       running_zscored(x, wins)   856   864   913    871   897   3198   100 a c        
##                                        running_sharpe(x, wins)   852   858   902    863   887   2939   100 a c        
##                                    running_apx_median(x, wins) 14707 14831 15502  14961 15805  20289   100         i  
##                                      running_centered(x, wins)   620   626   669    635   667   2466   100 a c        
##                                        running_scaled(x, wins)   851   858   886    869   893   1131   100 a c
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
|running_sum(x, wins)                                                         |       0.51|      1.88|        3.72|
|running_sd3(x, wins)                                                         |       8.72|     17.75|        2.03|
|sd3(x)                                                                       |      10.52|     18.91|        1.80|
|sd3(x, wts = w)                                                              |      13.06|     20.37|        1.56|
|as.centsums(x1, 2)                                                           |       1.83|      2.54|        1.38|
|running_tstat(x, wins)                                                       |       9.32|     11.50|        1.23|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.04|        1.19|
|running_scaled(x, wins)                                                      |       9.03|     10.75|        1.19|
|running_sharpe(x, wins)                                                      |       9.22|     10.94|        1.19|
|running_skew(x, wins)                                                        |      39.66|     46.81|        1.18|
|running_skew4(x, wins)                                                       |      48.07|     56.51|        1.18|
|running_zscored(x, wins)                                                     |       9.49|     11.08|        1.17|
|kurt5(x)                                                                     |     173.12|    200.94|        1.16|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.87|        1.16|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.04|        1.14|
|skew4(x)                                                                     |      97.66|    111.47|        1.14|
|mobj3 %-% mobj1                                                              |       0.26|      0.30|        1.14|
|running_kurt5(x, wins)                                                       |      76.25|     86.35|        1.13|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      6.83|        1.08|
|skewness(x)                                                                  |     100.16|    108.62|        1.08|
|kurtosis(x)                                                                  |      98.46|    105.03|        1.07|
|running_centered(x, wins)                                                    |       7.68|      8.12|        1.06|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      5.49|        1.06|
|running_apx_median(x, wins)                                                  |     179.28|    188.12|        1.05|
|dumbk(x)                                                                     |     209.85|    217.76|        1.04|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    145.70|        1.04|
|sd(x)                                                                        |       5.82|      5.95|        1.02|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    202.47|        1.01|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|mean(x)                                                                      |       2.07|      2.07|        1.00|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.99|
|running_kurt(x, wins)                                                        |      70.61|     67.58|        0.96|
|obj3 %-% obj1                                                                |       0.19|      0.18|        0.95|
|as.centsums(x1, 4)                                                           |      20.20|     19.25|        0.95|
|as.centsums(x1, 3)                                                           |      11.79|     11.21|        0.95|
|roll::roll_sd(xm, wins)                                                      |      77.75|     73.07|        0.94|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.47|        0.93|
|roll::roll_mean(xm, wins)                                                    |      28.39|     26.24|        0.92|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     23.86|        0.91|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     24.41|        0.91|
|roll::roll_sum(xm, wins)                                                     |      27.28|     24.55|        0.90|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    834.48|        0.88|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      2.78|        0.88|
|c(obj1, obj2)                                                                |       0.25|      0.22|        0.88|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    439.34|        0.88|
|as.centsums(x1, 1)                                                           |       2.83|      2.40|        0.85|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      3.62|        0.81|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      2.53|        0.68|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.40|        0.59|


