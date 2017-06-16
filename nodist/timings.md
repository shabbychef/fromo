

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
##                                "24ddfa411137"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.5010      
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
##       sum(x)    80    80    83     80    81   120   100 a      
##      mean(x)   162   164   174    165   167   282   100 a      
##        sd(x)   459   464   492    472   496   838   100  b     
##  skewness(x)  8269  8311  8723   8373  8813 10644   100    d   
##  kurtosis(x)  8096  8146  8486   8200  8533 11247   100    d   
##       sd3(x)  1463  1469  1535   1482  1549  2660   100   c    
##     skew4(x)  8835  8869  9182   8930  9192 12982   100     e  
##     kurt5(x) 15776 15852 16492  15938 16531 20232   100      f 
##     dumbk(x) 17033 17103 17993  17340 18720 23092   100       g
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16.2 16.3 17.1   16.9 17.6 20.0   100   b
##                                                               sd3(x, wts = w)  1.6  1.7  1.8    1.8  1.8  2.1   100  a 
##                                                                 slow_sd(x, w)  1.4  1.5  2.0    1.5  2.8  3.7   100  a
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
##  as.centsums(x1, 1)  195  199  204    201  206  280   100 a  
##  as.centsums(x1, 2)  199  201  208    203  206  426   100 a  
##  as.centsums(x1, 3)  900  903  921    905  926 1110   100  b 
##  as.centsums(x1, 4) 1541 1544 1579   1548 1570 2217   100   c
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
##  c(obj1, obj2)  14 15   18     16 16 179   100   b
##  obj3 %-% obj1  11 12   14     12 13  91   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  2.7    2.5 2.6 24.5   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.1  2.2    2.2 2.3  2.9   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.1  2.7    2.3 2.4 15.3   100   a
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
##  as.centcosums(x1, max_ord)  61 62   69     64 75 155   100   b
##             mobj3 %-% mobj1  19 21   24     22 24  38   100  a
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
##                                                           expr   min    lq  mean median    uq    max neval          cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32222 32854 35743  34547 35192 109642   100           k 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 66208 68649 71314  69641 70966 148701   100            l
##                                           running_sum(x, wins)   103   108   135    110   120   2069   100 a           
##                                          running_mean(x, wins)   105   109   115    112   119    150   100 a           
##                                       roll::roll_sum(xm, wins)  1831  1912  1992   1950  2007   3945   100   cd        
##                                      roll::roll_mean(xm, wins)  2015  2083  2168   2117  2186   3542   100    de       
##                                        roll::roll_sd(xm, wins)  5715  5757  6199   5812  5890  32573   100       gh    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   360   379   553    401   445   2702   100 a  d        
##                             RollingWindow::RollingSum(x, wins)   113   123   207    138   183   2175   100 a           
##                            RollingWindow::RollingMean(x, wins)   138   154   229    177   218   2352   100 ab          
##                             RollingWindow::RollingStd(x, wins)   225   243   315    259   316   2514   100 a c         
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1733  1807  1959   1949  1998   4003   100  bcd        
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1720  1829  2003   2003  2051   3780   100   cd        
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10090 10589 12170  10932 12605  85323   100         i   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   436   442   459    451   474    540   100 a  d        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   570   588   609    596   614    860   100 a  d        
##                                           running_sd3(x, wins)  1396  1444  1539   1471  1532   3589   100 a  d        
##                                          running_skew(x, wins)  3721  3747  3830   3763  3805   6024   100     ef      
##                                         running_skew4(x, wins)  4410  4473  4621   4536  4638   6633   100      fg     
##                                          running_kurt(x, wins)  5456  5483  5635   5521  5664   6943   100       gh    
##                                         running_kurt5(x, wins)  6777  6839  7099   6903  7012   9394   100        h    
##                                         running_tstat(x, wins)   891   900   944    904   939   2657   100 a  d        
##                                       running_zscored(x, wins)   858   866   938    876   912   3218   100 a  d        
##                                        running_sharpe(x, wins)   859   864   915    881   911   2847   100 a  d        
##                                    running_apx_median(x, wins) 14761 14865 15491  15022 16041  18140   100          j  
##                                      running_centered(x, wins)   629   637   685    647   667   2521   100 a  d        
##                                        running_scaled(x, wins)   856   862   887    869   892   1063   100 a  d
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
|running_sum(x, wins)                                                         |       0.51|      1.63|        3.23|
|running_sd3(x, wins)                                                         |       8.72|     18.65|        2.14|
|sd3(x)                                                                       |      10.52|     18.61|        1.77|
|sd3(x, wts = w)                                                              |      13.06|     21.30|        1.63|
|as.centsums(x1, 2)                                                           |       1.83|      2.52|        1.38|
|running_tstat(x, wins)                                                       |       9.32|     11.44|        1.23|
|running_sharpe(x, wins)                                                      |       9.22|     11.09|        1.20|
|running_zscored(x, wins)                                                     |       9.49|     11.37|        1.20|
|running_scaled(x, wins)                                                      |       9.03|     10.75|        1.19|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |       6.29|      7.39|        1.17|
|running_skew(x, wins)                                                        |      39.66|     46.42|        1.17|
|running_skew4(x, wins)                                                       |      48.07|     56.01|        1.17|
|kurt5(x)                                                                     |     173.12|    199.90|        1.15|
|skew4(x)                                                                     |      97.66|    111.30|        1.14|
|running_kurt5(x, wins)                                                       |      76.25|     86.04|        1.13|
|mobj3 %-% mobj1                                                              |       0.26|      0.29|        1.12|
|as.centcosums(x1, max_ord)                                                   |       0.75|      0.84|        1.12|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|        1.10|
|running_centered(x, wins)                                                    |       7.68|      8.30|        1.08|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |       5.19|      5.56|        1.07|
|skewness(x)                                                                  |     100.16|    105.73|        1.06|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |     140.59|    147.51|        1.05|
|running_apx_median(x, wins)                                                  |     179.28|    187.76|        1.05|
|kurtosis(x)                                                                  |      98.46|    102.86|        1.04|
|dumbk(x)                                                                     |     209.85|    218.09|        1.04|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    207.81|        1.04|
|sd(x)                                                                        |       5.82|      5.96|        1.02|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.03|        1.02|
|mean(x)                                                                      |       2.07|      2.11|        1.02|
|sum(x)                                                                       |       1.00|      1.00|        1.00|
|running_kurt(x, wins)                                                        |      70.61|     68.30|        0.97|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |       6.93|      6.70|        0.97|
|roll::roll_sd(xm, wins)                                                      |      77.75|     75.13|        0.97|
|as.centsums(x1, 4)                                                           |      20.20|     19.14|        0.95|
|as.centsums(x1, 3)                                                           |      11.79|     11.16|        0.95|
|roll::roll_mean(xm, wins)                                                    |      28.39|     26.28|        0.93|
|obj3 %-% obj1                                                                |       0.19|      0.17|        0.92|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |     943.85|    864.38|        0.92|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |      26.10|     23.75|        0.91|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |      26.80|     24.28|        0.91|
|c(obj1, obj2)                                                                |       0.25|      0.22|        0.89|
|roll::roll_sum(xm, wins)                                                     |      27.28|     24.15|        0.89|
|RollingWindow::RollingMean(x, wins)                                          |       3.15|      2.77|        0.88|
|as.centsums(x1, 1)                                                           |       2.83|      2.47|        0.87|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |     501.06|    433.24|        0.86|
|RollingWindow::RollingStd(x, wins)                                           |       4.48|      3.82|        0.85|
|unjoin_cent_sums(rs3, rs2)                                                   |       0.03|      0.03|        0.81|
|RollingWindow::RollingSum(x, wins)                                           |       3.72|      2.51|        0.67|
|running_sum(x, wins, robust = FALSE)                                         |       0.63|      0.42|        0.66|
|running_mean(x, wins)                                                        |       2.36|      1.40|        0.59|


