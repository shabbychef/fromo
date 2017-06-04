

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
##                                "a7b6c71d43d2"                                      "x86_64"                                     "unknown" 
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
##       sum(x)    80    82   101     93   111   290   100 a     
##      mean(x)   162   176   213    193   237   405   100 a     
##        sd(x)   459   512   578    557   634   916   100 ab    
##  skewness(x)  8376  9468 10206   9954 10739 15920   100   cd  
##  kurtosis(x)  8679  9862 10623  10587 11277 14299   100    d  
##       sd3(x)   856   941  1031   1002  1111  1359   100  b    
##     skew4(x)  8408  9417 10055  10006 10816 12291   100   c   
##     kurt5(x) 15012 16763 17939  17720 18975 24105   100     e 
##     dumbk(x) 17564 19787 21428  21238 22990 29235   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16304 17475 18270  18132 18655 27719   100   c
##                                                               sd3(x, wts = w)   981  1015  1111   1085  1163  1489   100 a  
##                                                                 slow_sd(x, w)  1398  1503  2120   1833  2661  4213   100  b
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
##  as.centsums(x1, 1)  188  194  215    200  219  387   100  b  
##  as.centsums(x1, 2)  113  118  131    122  133  422   100 a   
##  as.centsums(x1, 3)  856  900  957    921  960 1358   100   c 
##  as.centsums(x1, 4) 1481 1567 1662   1612 1686 2332   100    d
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
##  c(obj1, obj2)  14 16   23     19 25 196   100   a
##  obj3 %-% obj1  11 12   19     15 20 136   100   a
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
##    join_cent_sums(rs1, rs2) 2.1 2.2  2.6    2.3 2.5  25   100   a
##  unjoin_cent_sums(rs3, rs2) 1.8 2.0  2.3    2.1 2.3  16   100   a
##  unjoin_cent_sums(rs3, rs1) 1.8 2.0  2.4    2.1 2.3  14   100   a
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
##  as.centcosums(x1, max_ord)  54 56   63     57 69 162   100   b
##             mobj3 %-% mobj1  18 19   21     20 21  45   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32414 37694 40289  39440 41691 111644   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 65776 74544 80322  78869 81033 172658   100           k
##                                           running_sum(x, wins)    74    83   140     94   107   2724   100 a          
##                                          running_mean(x, wins)    73    87   121     95   111   2119   100 a          
##                                       roll::roll_sum(xm, wins)  1905  2022  2280   2189  2410   3842   100  b d       
##                                      roll::roll_mean(xm, wins)  2076  2246  2526   2465  2733   4165   100   cde      
##                                        roll::roll_sd(xm, wins)  5700  6044  6686   6587  7034   9721   100       g    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   359   424   572    478   551   3502   100 abc        
##                             RollingWindow::RollingSum(x, wins)   112   168   260    197   238   3379   100 ab         
##                            RollingWindow::RollingMean(x, wins)   140   180   323    229   281   3364   100 ab         
##                             RollingWindow::RollingStd(x, wins)   225   300   502    342   413   3263   100 abc        
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1769  2043  2254   2226  2385   3040   100  b d       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1739  2089  2291   2220  2507   3067   100  b d       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10366 12092 13015  12845 14026  16468   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   380   416   463    448   503    703   100 ab         
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   445   487   544    525   590    872   100 abc        
##                                           running_sd3(x, wins)   612   691  1652    751   829  83261   100 abc        
##                                          running_skew(x, wins)  3597  4011  4260   4202  4576   5249   100    d f     
##                                         running_skew4(x, wins)  3653  4002  4358   4249  4622   8099   100     ef     
##                                          running_kurt(x, wins)  5245  5683  6112   6077  6435   7436   100      fg    
##                                         running_kurt5(x, wins)  5875  6464  6939   6895  7369  11298   100       g    
##                                         running_tstat(x, wins)   669   727   826    790   875   2571   100 abc        
##                                       running_zscored(x, wins)   663   727   793    756   835   1174   100 abc        
##                                        running_sharpe(x, wins)   661   730   791    760   832   1127   100 abc        
##                                    running_apx_median(x, wins) 14088 15727 16760  16856 17431  21697   100         i  
##                                      running_centered(x, wins)   570   629   690    662   726    971   100 abc        
##                                        running_scaled(x, wins)   668   731   845    793   905   2901   100 abc
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
|running_sum(x, wins)                                                         |       0.51|      1.45|        2.86|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      4.21|        1.34|
|running_skew(x, wins)                                                        |      39.66|     50.08|        1.26|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      5.43|        1.21|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    167.33|        1.19|
|running_apx_median(x, wins)                                                  |     179.28|    199.89|        1.11|
|running_scaled(x, wins)                                                      |       9.03|      9.97|        1.10|
|running_sharpe(x, wins)                                                      |       9.22|     10.11|        1.10|
|running_centered(x, wins)                                                    |       7.68|      8.35|        1.09|
|running_kurt5(x, wins)                                                       |      76.25|     82.82|        1.09|
|running_skew4(x, wins)                                                       |      48.07|     51.87|        1.08|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    216.41|        1.08|
|skew4(x)                                                                     |      97.66|    104.96|        1.07|
|kurt5(x)                                                                     |     173.12|    185.38|        1.07|
|running_zscored(x, wins)                                                     |       9.49|     10.09|        1.06|
|skewness(x)                                                                  |     100.16|    106.39|        1.06|
|running_kurt(x, wins)                                                        |      70.61|     74.89|        1.06|
|running_tstat(x, wins)                                                       |       9.32|      9.84|        1.06|
|kurtosis(x)                                                                  |      98.46|    103.74|        1.05|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      5.46|        1.05|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     28.05|        1.05|
|dumbk(x)                                                                     |     209.85|    218.42|        1.04|
|sd(x)                                                                        |       5.82|      6.04|        1.04|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      6.51|        1.03|
|running_sd3(x, wins)                                                         |       8.72|      9.00|        1.03|
|sd3(x)                                                                       |      10.52|     10.78|        1.02|
|roll::roll_sd(xm, wins)                                                      |      77.75|     79.63|        1.02|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     26.71|        1.02|
|sd3(x, wts = w)                                                              |      13.06|     13.27|        1.02|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    956.41|        1.01|
|roll::roll_mean(xm, wins)                                                    |      28.39|     28.69|        1.01|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      3.75|        1.01|
|mean(x)                                                                      |       2.07|      2.08|        1.01|
|as.centsums(x1, 4)                                                           |      20.20|     20.28|        1.00|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|roll::roll_sum(xm, wins)                                                     |      27.28|     27.14|        0.99|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.82|        0.98|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    491.96|        0.98|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.73|        0.98|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        0.97|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        0.97|
|mobj3 %-% mobj1                                                              |       0.26|      0.25|        0.97|
|as.centsums(x1, 3)                                                           |      11.79|     11.24|        0.95|
|as.centsums(x1, 1)                                                           |       2.83|      2.66|        0.94|
|obj3 %-% obj1                                                                |       0.19|      0.17|        0.90|
|as.centsums(x1, 2)                                                           |       1.83|      1.63|        0.89|
|c(obj1, obj2)                                                                |       0.25|      0.22|        0.88|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.79|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.44|        0.61|


