

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
##                                "b8fa3e18d527"                                      "x86_64"                                     "unknown" 
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
##       sum(x)    79    80    82     80    81   134   100 a      
##      mean(x)   162   163   169    165   168   267   100 a      
##        sd(x)   459   465   489    477   491   889   100  b     
##  skewness(x)  8241  8283  8615   8394  8668 10788   100    de  
##  kurtosis(x)  8077  8128  8421   8212  8406 10111   100    d   
##       sd3(x)   849   859   884    863   876  1251   100   c    
##     skew4(x)  8335  8373  8683   8454  8653 20863   100     e  
##     kurt5(x) 14990 15053 15333  15161 15443 17015   100      f 
##     dumbk(x) 17048 17229 17977  17734 18626 20596   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16110 16225 16613  16431 16658 19580   100   c
##                                                               sd3(x, wts = w)   981   987  1021   1009  1026  1188   100 a  
##                                                                 slow_sd(x, w)  1389  1415  1933   1466  2774  3695   100  b
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
##  as.centsums(x1, 1)  188  189  199    191  202  260   100  b  
##  as.centsums(x1, 2)  114  115  122    117  128  287   100 a   
##  as.centsums(x1, 3)  851  855  896    871  905 1185   100   c 
##  as.centsums(x1, 4) 1483 1487 1557   1509 1594 2081   100    d
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
##  c(obj1, obj2)  14 15   18     15 16 191   100   a
##  obj3 %-% obj1  11 12   14     12 13  92   100   a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  2.7    2.4 2.6 28.1   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.2    2.1 2.2  3.1   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.7    2.1 2.3 29.0   100   a
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
##  as.centcosums(x1, max_ord)  51 54   60     55 65 149   100   b
##             mobj3 %-% mobj1  16 18   21     19 20  42   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32673 33920 37195  35430 38037 129476   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64588 66966 72316  68303 73202 142586   100           k
##                                           running_sum(x, wins)    72    77    87     81    95    160   100 a          
##                                          running_mean(x, wins)    73    78    90     84    96    183   100 a          
##                                       roll::roll_sum(xm, wins)  1785  1859  2016   1928  2070   3686   100 abc        
##                                      roll::roll_mean(xm, wins)  1961  2014  2187   2115  2307   3199   100  b d       
##                                        roll::roll_sd(xm, wins)  5436  5587  6073   5775  6249  11006   100      fg    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   362   394   547    432   485   3688   100 ab         
##                             RollingWindow::RollingSum(x, wins)   109   136   292    168   208   2308   100 ab         
##                            RollingWindow::RollingMean(x, wins)   138   165   985    209   242  77430   100 ab         
##                             RollingWindow::RollingStd(x, wins)   224   245   378    296   334   2560   100 ab         
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1763  1945  2135   2010  2176   4203   100 a  d       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1744  2003  2160   2071  2208   4211   100  b d       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10310 10929 12770  11431 12841  88209   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   379   390   450    415   453   2892   100 ab         
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   440   458   499    472   516    796   100 ab         
##                                           running_sd3(x, wins)   612   631   735    654   718   3804   100 ab         
##                                          running_skew(x, wins)  3570  3643  3842   3763  3907   4706   100   cde      
##                                         running_skew4(x, wins)  3657  3746  4096   3871  4296   6009   100    def     
##                                          running_kurt(x, wins)  5229  5382  5694   5527  5906   7426   100     e g    
##                                         running_kurt5(x, wins)  5907  6068  6510   6200  6800   9232   100       g    
##                                         running_tstat(x, wins)   665   675   788    699   770   3607   100 ab         
##                                       running_zscored(x, wins)   669   688  1153    719   789  41341   100 ab         
##                                        running_sharpe(x, wins)   661   674  1093    706   772  33931   100 ab         
##                                    running_apx_median(x, wins) 14103 14588 15524  15038 16202  20390   100         i  
##                                      running_centered(x, wins)   566   582   671    608   661   2808   100 ab         
##                                        running_scaled(x, wins)   660   674   737    701   766   1036   100 ab
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
allt %>% select(-sernum) %>% distinct(expr, .keep_all = TRUE) %>% 
    arrange(desc(last_status)) %>% head(n = 50) %>% 
    kable()
```



|expr                                                                         | meantime| sumx_time| normalized| first_norm| last_norm| relchange| last_status|
|:----------------------------------------------------------------------------|--------:|---------:|----------:|----------:|---------:|---------:|-----------:|
|running_sum(x, wins)                                                         |  4.3e+04|     48318|       0.88|       0.88|      1.45|         1|        1.63|
|running_skew(x, wins)                                                        |  3.3e+06|     83136|      39.66|      39.66|     50.08|         1|        1.26|
|RollingWindow::RollingStd(x, wins)                                           |  3.7e+05|     82037|       4.48|       4.48|      5.43|         1|        1.21|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |  1.2e+07|     82037|     140.59|     140.59|    167.33|         1|        1.19|
|running_kurt(x, wins)                                                        |  5.9e+06|     83136|      70.61|      70.61|     74.89|         1|        1.06|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |  4.3e+05|     82037|       5.19|       5.19|      5.46|         1|        1.05|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |  5.2e+05|     82037|       6.29|       6.29|      6.51|         1|        1.03|
|running_sd3(x, wins)                                                         |  7.2e+05|     82037|       8.72|       8.72|      9.00|         1|        1.03|
|sum(x)                                                                       |  4.8e+04|     48318|       1.00|       1.00|      1.00|         1|        1.00|
|RollingWindow::RollingMean(x, wins)                                          |  2.7e+05|     48318|       5.51|       5.51|      4.21|         1|        0.76|
|running_sum(x, wins, robust = FALSE)                                         |  5.3e+04|     85049|       0.63|       0.63|      0.42|         1|        0.66|
|running_apx_median(x, wins)                                                  |  1.5e+07|     48318|     313.88|     313.88|    199.89|         1|        0.64|
|running_scaled(x, wins)                                                      |  7.6e+05|     48318|      15.82|      15.82|      9.97|         1|        0.63|
|running_sharpe(x, wins)                                                      |  7.8e+05|     48318|      16.14|      16.14|     10.11|         1|        0.63|
|running_centered(x, wins)                                                    |  6.5e+05|     48318|      13.44|      13.44|      8.35|         1|        0.62|
|running_kurt5(x, wins)                                                       |  6.5e+06|     48318|     133.50|     133.50|     82.82|         1|        0.62|
|running_skew4(x, wins)                                                       |  4.1e+06|     48318|      84.16|      84.16|     51.87|         1|        0.62|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |  1.7e+07|     48318|     351.31|     351.31|    216.41|         1|        0.62|
|skew4(x)                                                                     |  8.3e+06|     48318|     170.98|     170.98|    104.96|         1|        0.61|
|kurt5(x)                                                                     |  1.5e+07|     48318|     303.10|     303.10|    185.38|         1|        0.61|
|running_zscored(x, wins)                                                     |  8.0e+05|     48318|      16.62|      16.62|     10.09|         1|        0.61|
|skewness(x)                                                                  |  8.5e+06|     48318|     175.37|     175.37|    106.39|         1|        0.61|
|running_tstat(x, wins)                                                       |  7.9e+05|     48318|      16.31|      16.31|      9.84|         1|        0.60|
|kurtosis(x)                                                                  |  8.3e+06|     48318|     172.39|     172.39|    103.74|         1|        0.60|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |  2.3e+06|     48318|      46.93|      46.93|     28.05|         1|        0.60|
|dumbk(x)                                                                     |  1.8e+07|     48318|     367.42|     367.42|    218.42|         1|        0.59|
|sd(x)                                                                        |  4.9e+05|     48318|      10.19|      10.19|      6.04|         1|        0.59|
|sd3(x)                                                                       |  8.9e+05|     48318|      18.42|      18.42|     10.78|         1|        0.59|
|roll::roll_sd(xm, wins)                                                      |  6.6e+06|     48318|     136.13|     136.13|     79.63|         1|        0.59|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |  2.2e+06|     48318|      45.69|      45.69|     26.71|         1|        0.58|
|sd3(x, wts = w)                                                              |  1.1e+06|     48318|      22.86|      22.86|     13.27|         1|        0.58|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |  8.0e+07|     48318|    1652.53|    1652.53|    956.41|         1|        0.58|
|roll::roll_mean(xm, wins)                                                    |  2.4e+06|     48318|      49.71|      49.71|     28.69|         1|        0.58|
|RollingWindow::RollingSum(x, wins)                                           |  3.1e+05|     48318|       6.51|       6.51|      3.75|         1|        0.58|
|mean(x)                                                                      |  1.7e+05|     48318|       3.62|       3.62|      2.08|         1|        0.57|
|as.centsums(x1, 4)                                                           |  1.7e+06|     48318|      35.36|      35.36|     20.28|         1|        0.57|
|roll::roll_sum(xm, wins)                                                     |  2.3e+06|     48318|      47.77|      47.77|     27.14|         1|        0.57|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |  5.9e+05|     48318|      12.13|      12.13|      6.82|         1|        0.56|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |  4.2e+07|     48318|     877.27|     877.27|    491.96|         1|        0.56|
|as.centcosums(x1, max_ord)                                                   |  6.4e+04|     48318|       1.32|       1.32|      0.73|         1|        0.56|
|join_cent_sums(rs1, rs2)                                                     |  2.7e+03|     48318|       0.06|       0.06|      0.03|         1|        0.56|
|unjoin_cent_sums(rs3, rs1)                                                   |  2.5e+03|     48318|       0.05|       0.05|      0.03|         1|        0.55|
|mobj3 %-% mobj1                                                              |  2.2e+04|     48318|       0.45|       0.45|      0.25|         1|        0.55|
|as.centsums(x1, 3)                                                           |  1.0e+06|     48318|      20.65|      20.65|     11.24|         1|        0.54|
|as.centsums(x1, 1)                                                           |  2.4e+05|     48318|       4.96|       4.96|      2.66|         1|        0.54|
|obj3 %-% obj1                                                                |  1.6e+04|     48318|       0.33|       0.33|      0.17|         1|        0.51|
|as.centsums(x1, 2)                                                           |  1.6e+05|     48318|       3.21|       3.21|      1.63|         1|        0.51|
|c(obj1, obj2)                                                                |  2.1e+04|     48318|       0.43|       0.43|      0.22|         1|        0.50|
|running_sum(x, wins, robust = TRUE)                                          |  1.6e+05|     85049|       1.88|       1.88|      0.90|         1|        0.48|
|unjoin_cent_sums(rs3, rs2)                                                   |  2.8e+03|     48318|       0.06|       0.06|      0.03|         1|        0.45|


