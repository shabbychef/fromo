

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
##                                "261486475bd6"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4400      
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
##         expr   min    lq  mean median    uq   max neval   cld
##       sum(x)    80    80    97     87    99   287   100 a    
##      mean(x)   162   166   195    179   211   343   100 a    
##        sd(x)   460   469   544    504   570  1024   100 ab   
##  skewness(x)  8301  8665  9990   9923 10865 15263   100   c  
##  kurtosis(x)  8120  8237  9614   9627 10459 13230   100   c  
##       sd3(x)   848   859   986    932  1097  1328   100  b   
##     skew4(x)  8575  8653  9840   9397 10754 15992   100   c  
##     kurt5(x) 15403 15829 17920  18068 19515 22010   100    d 
##     dumbk(x) 17055 18634 20368  19961 22205 29449   100     e
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16297 16451 17165  16620 17377 24266   100   c
##                                                               sd3(x, wts = w)   980   984  1024   1001  1018  1741   100 a  
##                                                                 slow_sd(x, w)  1391  1421  1952   1518  2597  3708   100  b
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
##  as.centsums(x1, 1)  193  209  229    211  237  412   100  b  
##  as.centsums(x1, 2)  114  124  144    126  144  321   100 a   
##  as.centsums(x1, 3)  862  936  990    941 1010 1349   100   c 
##  as.centsums(x1, 4) 1503 1636 1785   1681 1959 2346   100    d
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
##  c(obj1, obj2)  15 18   26     21 29 257   100   a
##  obj3 %-% obj1  12 14   23     17 26 153   100   a
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
##    join_cent_sums(rs1, rs2) 2.5 2.7  3.9    2.9 3.2  51   100   a
##  unjoin_cent_sums(rs3, rs2) 2.3 2.5  3.9    2.7 2.9  75   100   a
##  unjoin_cent_sums(rs3, rs1) 2.3 2.5  3.3    2.7 3.0  47   100   a
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
##  as.centcosums(x1, max_ord)  56 59   73     61 79 187   100   b
##             mobj3 %-% mobj1  18 20   25     21 24  79   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 31793 33275 36183  34947 36669 108920   100         i 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 65520 68439 74873  70756 75372 145223   100          j
##                                           running_sum(x, wins)    73    78    88     82    94    137   100 a         
##                                          running_mean(x, wins)    73    77    88     84    94    140   100 a         
##                                       roll::roll_sum(xm, wins)  1890  1946  2049   2000  2077   2873   100  bcd      
##                                      roll::roll_mean(xm, wins)  2074  2112  2263   2179  2356   3022   100    de     
##                                        roll::roll_sd(xm, wins)  5737  5800  6179   5894  6479   8205   100      f    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   355   381   518    412   471   2597   100 a  d      
##                             RollingWindow::RollingSum(x, wins)   110   134   200    155   192   2020   100 a         
##                            RollingWindow::RollingMean(x, wins)   141   161   345    187   238   2811   100 ab        
##                             RollingWindow::RollingStd(x, wins)   224   253   400    284   324   2576   100 a c       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1756  1943  2040   1988  2119   2752   100  bcd      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1723  1993  2124   2045  2194   4518   100   cd      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10184 10728 11780  11456 12548  16826   100       g   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   355   365   396    383   398    651   100 a c       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   413   419   457    430   461    724   100 a c       
##                                           running_sd3(x, wins)   584   595  1427    616   666  73148   100 a  d      
##                                          running_skew(x, wins)  3630  3670  3892   3766  3995   5639   100     e     
##                                         running_skew4(x, wins)  3644  3683  3924   3770  4010   5747   100     e     
##                                          running_kurt(x, wins)  5328  5391  5712   5575  5884   7425   100      f    
##                                         running_kurt5(x, wins)  5942  6002  6419   6138  6633   8293   100      f    
##                                         running_tstat(x, wins)   622   632   687    654   700   1201   100 a  d      
##                                       running_zscored(x, wins)   625   630   706    652   698   2881   100 a  d      
##                                        running_sharpe(x, wins)   620   632   687    663   701    897   100 a  d      
##                                    running_apx_median(x, wins) 14065 14321 15247  14689 16232  18678   100        h  
##                                      running_centered(x, wins)   556   567   599    583   622    736   100 a  d      
##                                        running_scaled(x, wins)   620   631   698    656   701   2404   100 a  d
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
|RollingWindow::RollingSum(x, wins)                                           |       3.72|     12.49|        3.36|
|running_sum(x, wins)                                                         |       0.51|      0.98|        1.95|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      3.83|        1.22|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.93|        1.10|
|running_skew(x, wins)                                                        |      39.66|     43.13|        1.09|
|kurt5(x)                                                                     |     173.12|    183.06|        1.06|
|skew4(x)                                                                     |      97.66|    102.73|        1.05|
|mobj3 %-% mobj1                                                              |       0.26|      0.27|        1.04|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.76|        1.02|
|mean(x)                                                                      |       2.07|      2.09|        1.01|
|skewness(x)                                                                  |     100.16|    100.82|        1.01|
|kurtosis(x)                                                                  |      98.46|     98.79|        1.00|
|sd3(x)                                                                       |      10.52|     10.53|        1.00|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|dumbk(x)                                                                     |     209.85|    208.23|        0.99|
|sd(x)                                                                        |       5.82|      5.75|        0.99|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        0.99|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        0.97|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.68|        0.96|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    191.57|        0.95|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    134.07|        0.95|
|running_apx_median(x, wins)                                                  |     179.28|    166.76|        0.93|
|running_kurt5(x, wins)                                                       |      76.25|     70.88|        0.93|
|running_skew4(x, wins)                                                       |      48.07|     44.28|        0.92|
|running_centered(x, wins)                                                    |       7.68|      7.07|        0.92|
|running_kurt(x, wins)                                                        |      70.61|     63.91|        0.91|
|sd3(x, wts = w)                                                              |      13.06|     11.73|        0.90|
|as.centsums(x1, 4)                                                           |      20.20|     18.05|        0.89|
|running_scaled(x, wins)                                                      |       9.03|      7.99|        0.88|
|roll::roll_mean(xm, wins)                                                    |      28.39|     25.02|        0.88|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    828.47|        0.88|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      5.50|        0.87|
|obj3 %-% obj1                                                                |       0.19|      0.16|        0.87|
|as.centsums(x1, 3)                                                           |      11.79|     10.28|        0.87|
|c(obj1, obj2)                                                                |       0.25|      0.21|        0.87|
|running_sd3(x, wins)                                                         |       8.72|      7.50|        0.86|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     22.42|        0.86|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      4.45|        0.86|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     22.90|        0.85|
|roll::roll_sd(xm, wins)                                                      |      77.75|     66.40|        0.85|
|running_sharpe(x, wins)                                                      |       9.22|      7.76|        0.84|
|roll::roll_sum(xm, wins)                                                     |      27.28|     22.92|        0.84|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    412.92|        0.82|
|running_tstat(x, wins)                                                       |       9.32|      7.61|        0.82|
|as.centsums(x1, 1)                                                           |       2.83|      2.29|        0.81|
|running_zscored(x, wins)                                                     |       9.49|      7.59|        0.80|
|as.centsums(x1, 2)                                                           |       1.83|      1.43|        0.78|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.77|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.24|        0.53|


