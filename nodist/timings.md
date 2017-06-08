

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
##                                "850b5bb3499b"                                      "x86_64"                                     "unknown" 
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
##         expr   min    lq  mean median    uq   max neval     cld
##       sum(x)    80    80    81     80    81    90   100 a      
##      mean(x)   162   164   169    165   168   278   100 a      
##        sd(x)   459   462   475    468   481   602   100  b     
##  skewness(x)  8278  8313  8595   8374  8533 10489   100     e  
##  kurtosis(x)  8111  8131  8388   8187  8348 10466   100    d   
##       sd3(x)   848   851   870    860   881  1064   100   c    
##     skew4(x)  8535  8564  8654   8596  8706  9554   100     e  
##     kurt5(x) 15373 15426 15651  15493 15731 17211   100      f 
##     dumbk(x) 17062 17120 17712  17270 18133 20193   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16272 16348 16517  16437 16566 17503   100   c
##                                                               sd3(x, wts = w)   980   982  1005    991  1017  1215   100 a  
##                                                                 slow_sd(x, w)  1397  1409  1912   1449  2749  3730   100  b
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
##  as.centsums(x1, 1)  193  194  199    195  197  240   100  b  
##  as.centsums(x1, 2)  113  115  120    116  119  272   100 a   
##  as.centsums(x1, 3)  863  868  878    870  882  997   100   c 
##  as.centsums(x1, 4) 1506 1509 1524   1513 1533 1594   100    d
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
##  c(obj1, obj2)  14 15   19     16 17 179   100   a
##  obj3 %-% obj1  11 12   15     12 13 114   100   a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  2.9    2.3 2.5 26.3   100   b
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.1    2.1 2.2  2.9   100  a 
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.3    2.1 2.3 14.5   100  ab
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
##  as.centcosums(x1, max_ord)  53 54   60     56 62 145   100   b
##             mobj3 %-% mobj1  17 18   21     19 20  36   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32340 33174 35256  34692 35489 107522   100         i 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64760 66221 70167  67948 69411 144029   100          j
##                                           running_sum(x, wins)    73    77    83     79    86    120   100 a         
##                                          running_mean(x, wins)    74    77   107     80    88   2364   100 a         
##                                       roll::roll_sum(xm, wins)  1782  1827  2703   1869  1945  80413   100    de     
##                                      roll::roll_mean(xm, wins)  1965  2003  2092   2039  2129   3647   100   c e     
##                                        roll::roll_sd(xm, wins)  5427  5491  5638   5535  5623   7163   100      f    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   365   388   481    414   451   2600   100 a c       
##                             RollingWindow::RollingSum(x, wins)   112   138   255    151   192   2801   100 a         
##                            RollingWindow::RollingMean(x, wins)   144   163   313    183   222   2994   100 ab        
##                             RollingWindow::RollingStd(x, wins)   226   253   452    270   318   3301   100 a c       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1753  1945  2003   1980  2029   3968   100  bcd      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1731  1864  1989   2010  2055   2885   100  bcd      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10359 10639 11334  10818 12218  15091   100       g   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   356   361   378    370   387    479   100 a c       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   413   421   441    429   443    679   100 a c       
##                                           running_sd3(x, wins)   587   594   667    605   627   3061   100 a c       
##                                          running_skew(x, wins)  3635  3657  3725   3689  3759   4267   100     e     
##                                         running_skew4(x, wins)  3646  3675  3763   3704  3774   5886   100     e     
##                                          running_kurt(x, wins)  5329  5369  5518   5405  5563   6923   100      f    
##                                         running_kurt5(x, wins)  5937  5972  6143   6030  6139   8581   100      f    
##                                         running_tstat(x, wins)   621   630   651    641   653    880   100 a c       
##                                       running_zscored(x, wins)   624   630   680    647   668   2890   100 a c       
##                                        running_sharpe(x, wins)   621   629   672    644   660   2837   100 a c       
##                                    running_apx_median(x, wins) 14036 14131 14422  14238 14499  16871   100        h  
##                                      running_centered(x, wins)   556   565   586    575   591    747   100 a c       
##                                        running_scaled(x, wins)   618   628   671    643   659   2911   100 a c
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
|running_sum(x, wins)                                                         |       0.51|      0.96|        1.89|
|running_sd3(x, wins)                                                         |       8.72|     11.76|        1.35|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.04|        1.18|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      3.70|        1.17|
|obj3 %-% obj1                                                                |       0.19|      0.22|        1.15|
|kurt5(x)                                                                     |     173.12|    188.56|        1.09|
|running_skew(x, wins)                                                        |      39.66|     42.78|        1.08|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        1.07|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      4.79|        1.07|
|skew4(x)                                                                     |      97.66|    103.88|        1.06|
|skewness(x)                                                                  |     100.16|    104.39|        1.04|
|kurtosis(x)                                                                  |      98.46|    101.11|        1.03|
|c(obj1, obj2)                                                                |       0.25|      0.25|        1.02|
|dumbk(x)                                                                     |     209.85|    213.89|        1.02|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        1.01|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|mobj3 %-% mobj1                                                              |       0.26|      0.26|        1.00|
|sd3(x)                                                                       |      10.52|     10.43|        0.99|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.74|        0.99|
|mean(x)                                                                      |       2.07|      2.04|        0.99|
|sd(x)                                                                        |       5.82|      5.72|        0.98|
|roll::roll_sum(xm, wins)                                                     |      27.28|     26.71|        0.98|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    189.18|        0.94|
|running_apx_median(x, wins)                                                  |     179.28|    166.64|        0.93|
|running_kurt5(x, wins)                                                       |      76.25|     70.56|        0.93|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    129.83|        0.92|
|as.centsums(x1, 4)                                                           |      20.20|     18.58|        0.92|
|running_skew4(x, wins)                                                       |      48.07|     43.18|        0.90|
|running_kurt(x, wins)                                                        |      70.61|     63.08|        0.89|
|as.centsums(x1, 3)                                                           |      11.79|     10.49|        0.89|
|sd3(x, wts = w)                                                              |      13.06|     11.39|        0.87|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     22.71|        0.87|
|running_centered(x, wins)                                                    |       7.68|      6.66|        0.87|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    814.66|        0.86|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     23.10|        0.86|
|roll::roll_mean(xm, wins)                                                    |      28.39|     24.46|        0.86|
|roll::roll_sd(xm, wins)                                                      |      77.75|     66.37|        0.85|
|running_scaled(x, wins)                                                      |       9.03|      7.69|        0.85|
|as.centsums(x1, 1)                                                           |       2.83|      2.40|        0.85|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      4.35|        0.84|
|running_sharpe(x, wins)                                                      |       9.22|      7.63|        0.83|
|running_zscored(x, wins)                                                     |       9.49|      7.78|        0.82|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      5.61|        0.81|
|as.centsums(x1, 2)                                                           |       1.83|      1.48|        0.81|
|running_tstat(x, wins)                                                       |       9.32|      7.52|        0.81|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      5.04|        0.80|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    401.28|        0.80|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      2.55|        0.69|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_sum(x, wins, robust = TRUE)                                          |       1.88|      0.90|        0.48|


