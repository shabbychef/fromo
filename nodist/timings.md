

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
##                                               sysname                                               release 
##                                               "Linux"                                   "3.19.0-80-generic" 
##                                               version                                              nodename 
## "#88~14.04.1-Ubuntu SMP Fri Jan 13 14:54:07 UTC 2017"                                        "3d9d358bf8bc" 
##                                               machine                                                 login 
##                                              "x86_64"                                             "unknown" 
##                                                  user                                        effective_user 
##                                                "spav"                                                "spav"
```

```r
print(sessionInfo())
```

```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS: /usr/lib/libblas/libblas.so.3.7.0
## LAPACK: /usr/lib/lapack/liblapack.so.3.7.0
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
##  [7] dplyr_0.7.0            ggplot2_2.2.1          knitr_1.16             fromo_0.1.3.6200      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11     magrittr_1.5     grDevices_3.4.0  munsell_0.4.3    colorspace_1.3-2 R6_2.2.1         rlang_0.1.1      stringr_1.2.0    plyr_1.8.4      
## [10] tools_3.4.0      grid_3.4.0       gtable_0.2.0     stats_3.4.0      lazyeval_0.2.0   assertthat_0.2.0 tibble_1.3.3     formatR_1.5      graphics_3.4.0  
## [19] glue_1.1.0       evaluate_0.10    stringi_1.1.5    compiler_3.4.0   scales_0.4.1
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
##       sum(x)    76    76    77     76    77    81   100 a     
##      mean(x)   154   155   158    156   157   246   100 a     
##        sd(x)   294   299   310    304   314   395   100 a     
##  skewness(x)  7278  7308  7540   7330  7436 10086   100   c   
##  kurtosis(x)  7156  7207  7839   7246  7414 45309   100   cd  
##       sd3(x)  1058  1062  1080   1067  1077  1300   100  b    
##     skew4(x)  7994  8014  8129   8045  8104  9466   100    d  
##     kurt5(x) 14273 14299 14478  14366 14479 16241   100     e 
##     dumbk(x) 14904 15025 15425  15159 15506 17555   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 14606 14621 14688  14644 14698 15647   100   b
##                                                               sd3(x, wts = w)  1058  1061  1072   1065  1077  1222   100  a 
##                                                                 slow_sd(x, w)   887   903  2023    912  1729 73022   100  a
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
##  as.centsums(x1, 1)   73   74   75     75   75   80   100 a   
##  as.centsums(x1, 2)  134  135  137    135  136  281   100  b  
##  as.centsums(x1, 3)  813  815  819    816  819  965   100   c 
##  as.centsums(x1, 4) 1410 1412 1415   1413 1416 1461   100    d
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
##  c(obj1, obj2)  15 15   18     16 16 183   100   b
##  obj3 %-% obj1  11 12   14     13 13  89   100  a
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
##    join_cent_sums(rs1, rs2) 2.4 2.5  2.8    2.6 2.7 23.9   100   b
##  unjoin_cent_sums(rs3, rs2) 2.1 2.3  2.4    2.4 2.5  2.9   100  a 
##  unjoin_cent_sums(rs3, rs1) 2.2 2.3  2.5    2.4 2.5  6.0   100  ab
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
##  as.centcosums(x1, max_ord)  49 50   53     51 51 136   100   b
##             mobj3 %-% mobj1  16 17   19     18 19  31   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 22357 22793 25804  23015 24434 104468   100         i 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 51702 52508 54788  52913 54331 128530   100          j
##                                           running_sum(x, wins)    99   102   106    103   105    154   100 a         
##                                          running_mean(x, wins)    97   102   106    104   106    142   100 a         
##                                       roll::roll_sum(xm, wins)  1647  1669  1723   1699  1739   3495   100 a c       
##                                      roll::roll_mean(xm, wins)  1803  1824  1864   1848  1891   2121   100   cd      
##                                        roll::roll_sd(xm, wins)  4693  4721  4797   4764  4824   6931   100     ef    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   387   402   537    418   466   2200   100 a c       
##                             RollingWindow::RollingSum(x, wins)   104   116   157    128   164   1813   100 ab        
##                            RollingWindow::RollingMean(x, wins)   133   144   224    159   207   1946   100 a c       
##                             RollingWindow::RollingStd(x, wins)   215   228   300    241   290   1779   100 a c       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1641  1720  1830   1844  1880   3197   100  bcd      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1638  1876  1856   1894  1911   2023   100   cd      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)  9223  9755 10131   9947 10116  12044   100       g   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   370   376   399    377   380   2093   100 a c       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   416   424   430    426   431    464   100 a c       
##                                           running_sd3(x, wins)  1245  1260  1356   1270  1289   2941   100 a c       
##                                          running_skew(x, wins)  3398  3406  3427   3412  3439   3613   100    de     
##                                         running_skew4(x, wins)  4077  4100  4144   4116  4148   5662   100     e     
##                                          running_kurt(x, wins)  4927  4939  4976   4949  4969   6636   100     ef    
##                                         running_kurt5(x, wins)  6219  6255  6365   6280  6309   7986   100      f    
##                                         running_tstat(x, wins)   833   838   845    840   848    946   100 a c       
##                                       running_zscored(x, wins)   787   792   801    794   804    942   100 a c       
##                                        running_sharpe(x, wins)   783   787   795    789   794    850   100 a c       
##                                    running_apx_median(x, wins) 13690 13792 13977  13860 13946  15609   100        h  
##                                      running_centered(x, wins)   550   555   561    557   564    587   100 a c       
##                                        running_scaled(x, wins)   783   787   796    790   794   1094   100 a c
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
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + labs(x = "release", y = "mean time taken, relative to sum(x)", 
    title = "fromo microbenchmark timings")
print(ph)
```

<img src="figure/all_timing_stats-1.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="600px" height="500px" />

```r
ph <- allt %>% ggplot(aes(sernum, relchange, group = expr, 
    color = expr)) + geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + labs(x = "release", y = "normalized time taken, relative to first iteration", 
    title = "fromo microbenchmark timings")
print(ph)
```

<img src="figure/all_timing_stats-2.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="600px" height="500px" />

```r
ph <- allt %>% ggplot(aes(sernum, relchange)) + geom_boxplot(aes(group = sernum), 
    alpha = 0.7) + stat_summary(aes(group = "1", color = "mean"), 
    fun.y = mean, geom = "line") + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    labs(x = "release", y = "normalized time taken, relative to first iteration", 
        color = "stat", title = "fromo microbenchmark timings")
print(ph)
```

<img src="figure/all_timing_stats-3.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="600px" height="500px" />

```r
allt %>% select(expr, sernum, relchange, last_status) %>% 
    tidyr::spread(key = "sernum", value = "relchange") %>% 
    arrange(desc(last_status)) %>% select(-last_status) %>% 
    head(n = 40) %>% kable()
```



|expr                                                                         | 0.1.3.3001| 0.1.3.3200| 0.1.3.3330| 0.1.3.3400| 0.1.3.4000| 0.1.3.4100| 0.1.3.4200| 0.1.3.4300| 0.1.3.4400| 0.1.3.4500| 0.1.3.4510| 0.1.3.5000| 0.1.3.5010| 0.1.3.5020| 0.1.3.6000| 0.1.3.6100| 0.1.3.6200|
|:----------------------------------------------------------------------------|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
|running_sum(x, wins)                                                         |          1|         NA|         NA|         NA|       2.72|       2.18|       1.95|         NA|       2.02|       2.53|       2.74|       2.45|       3.23|       3.72|       2.50|       2.70|       2.76|
|running_sd3(x, wins)                                                         |         NA|         NA|         NA|         NA|       1.00|       1.04|       0.86|         NA|       0.95|       0.78|       0.88|       0.80|       2.14|       2.03|       1.88|       2.05|       2.07|
|sd3(x)                                                                       |          1|       1.03|       1.01|       1.01|       1.03|       1.03|       1.02|       1.00|       1.02|       1.01|       0.99|       1.00|       1.77|       1.80|       1.71|       1.76|       1.33|
|join_cent_sums(rs1, rs2)                                                     |          1|       1.12|       0.96|       0.96|       1.02|       1.25|       1.12|       0.97|       1.10|       1.03|       1.05|       0.93|       1.02|       1.14|       1.01|       1.11|       1.19|
|running_skew4(x, wins)                                                       |          1|       0.94|       0.89|       0.89|       1.02|       1.03|       0.92|         NA|       0.97|       0.87|       0.96|       0.86|       1.17|       1.18|       1.04|       1.17|       1.16|
|running_tstat(x, wins)                                                       |          1|       0.92|       0.94|       0.94|       0.99|       1.01|       0.82|         NA|       0.86|       0.80|       0.90|       0.76|       1.23|       1.23|       1.10|       1.14|       1.16|
|running_skew(x, wins)                                                        |         NA|         NA|         NA|         NA|         NA|       1.00|       1.09|         NA|       1.16|       1.05|       1.18|       1.04|       1.17|       1.18|       1.05|       1.18|       1.15|
|running_scaled(x, wins)                                                      |          1|       0.94|       0.91|       0.91|       0.98|       1.03|       0.88|         NA|       0.92|       0.79|       0.88|       0.79|       1.19|       1.19|       1.08|       1.19|       1.15|
|running_zscored(x, wins)                                                     |          1|       0.90|       0.88|       0.88|       1.01|       0.98|       0.80|         NA|       0.89|       0.77|       0.84|       0.78|       1.20|       1.17|       1.02|       1.14|       1.12|
|running_sharpe(x, wins)                                                      |          1|       0.99|       0.90|       0.90|       0.98|       1.02|       0.84|         NA|       0.90|       0.79|       0.88|       0.77|       1.20|       1.19|       1.05|       1.11|       1.12|
|kurt5(x)                                                                     |          1|       1.02|       1.00|       1.00|       1.09|       1.10|       1.07|       1.06|       1.12|       1.09|       1.08|       1.08|       1.15|       1.16|       1.12|       1.14|       1.12|
|running_kurt5(x, wins)                                                       |          1|       0.94|       0.91|       0.91|       1.01|       1.03|       0.93|         NA|       1.00|       0.90|       0.98|       0.89|       1.13|       1.13|       1.01|       1.13|       1.11|
|skew4(x)                                                                     |          1|       1.01|       1.01|       1.01|       1.07|       1.07|       1.06|       1.05|       1.10|       1.08|       1.07|       1.09|       1.14|       1.14|       1.12|       1.13|       1.11|
|unjoin_cent_sums(rs3, rs1)                                                   |          1|       1.21|       0.97|       0.97|       0.91|       1.34|       0.88|       0.99|       0.98|       0.99|       0.98|       0.85|       1.10|       1.19|       0.99|       1.11|       1.11|
|sd3(x, wts = w)                                                              |          1|       0.94|       0.94|       0.94|       0.96|       0.97|       0.90|       0.90|       0.95|       0.90|       0.97|       0.83|       1.63|       1.56|       1.50|       1.54|       1.08|
|kurtosis(x)                                                                  |          1|       1.00|       1.00|       1.00|       1.06|       1.06|       1.03|       1.00|       1.05|       1.04|       1.03|       1.03|       1.04|       1.07|       1.04|       1.04|       1.04|
|running_apx_median(x, wins)                                                  |          1|       0.92|       0.90|       0.90|       1.06|       1.09|       0.93|         NA|       1.00|       0.89|       1.00|       0.88|       1.05|       1.05|       0.93|       1.05|       1.03|
|skewness(x)                                                                  |          1|       1.00|       1.00|       1.00|       1.06|       1.07|       1.04|       1.01|       1.06|       1.05|       1.03|       1.03|       1.06|       1.08|       1.06|       1.05|       1.03|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |          1|       0.93|       0.94|       0.94|       1.01|       1.07|       0.96|       0.95|       1.02|       0.96|       1.06|       0.91|       1.04|       1.01|       0.97|       1.00|       1.01|
|dumbk(x)                                                                     |          1|       1.03|       1.01|       1.01|       1.04|       1.04|       1.01|       0.99|       1.04|       1.02|       1.02|       1.02|       1.04|       1.04|       1.03|       1.03|       1.00|
|mean(x)                                                                      |          1|       1.02|       1.00|       1.00|       1.01|       1.01|       1.01|       1.01|       1.01|       1.06|       0.97|       1.03|       1.02|       1.00|       1.05|       1.02|       1.00|
|sum(x)                                                                       |          1|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |         NA|         NA|         NA|         NA|       1.00|       1.03|       0.86|         NA|       0.90|       0.85|       0.92|       0.83|       1.07|       1.06|       1.07|       1.02|       0.99|
|sd(x)                                                                        |          1|       1.02|       1.01|       1.01|       1.01|       1.01|       1.03|       0.99|       1.01|       1.02|       1.00|       1.02|       1.02|       1.02|       1.04|       1.02|       0.98|
|c(obj1, obj2)                                                                |          1|       1.29|       0.97|       0.97|       0.93|       0.98|       1.04|       0.87|       0.95|       0.86|       0.85|       0.81|       0.89|       0.88|       0.83|       0.96|       0.97|
|running_centered(x, wins)                                                    |          1|       0.93|       0.91|       0.91|       1.01|       0.98|       0.92|         NA|       0.94|       0.84|       0.91|       0.87|       1.08|       1.06|       0.94|       1.04|       0.97|
|RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)                  |         NA|         NA|         NA|         NA|       1.00|       1.11|       0.95|         NA|       1.00|       0.91|       1.01|       0.88|       1.05|       1.04|       0.91|       1.02|       0.96|
|obj3 %-% obj1                                                                |          1|       1.53|       0.96|       0.96|       0.91|       0.97|       1.00|       0.87|       1.00|       0.85|       0.85|       0.79|       0.92|       0.95|       0.90|       0.90|       0.96|
|running_kurt(x, wins)                                                        |         NA|         NA|         NA|         NA|         NA|       1.00|       0.91|         NA|       0.97|       0.85|       0.95|       0.86|       0.97|       0.96|       0.86|       0.96|       0.95|
|as.centsums(x1, 2)                                                           |          1|       0.82|       0.81|       0.81|       0.81|       0.88|       0.80|       0.78|       0.81|       0.77|       0.78|       0.72|       1.38|       1.38|       1.41|       1.38|       0.94|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |         NA|         NA|         NA|         NA|       1.00|       1.03|       0.87|         NA|       0.87|       0.77|       0.87|       0.75|       1.17|       1.08|       1.07|       0.96|       0.93|
|RollingWindow::RollingStd(x, wins)                                           |         NA|         NA|         NA|         NA|       1.00|       1.00|       1.10|         NA|       1.25|       0.98|       1.01|       1.00|       0.85|       0.81|       0.91|       1.08|       0.92|
|RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)                |          1|       0.92|       0.91|       0.91|       0.99|       1.01|       0.85|         NA|       0.92|       0.82|       0.91|       0.80|       0.91|       0.91|       0.83|       0.89|       0.92|
|as.centcosums(x1, max_ord)                                                   |          1|       0.93|       1.02|       1.02|       0.96|       1.24|       0.98|       1.02|       0.98|       0.87|       1.13|       0.85|       1.12|       1.16|       0.83|       0.93|       0.92|
|mobj3 %-% mobj1                                                              |          1|       0.93|       1.03|       1.03|       0.97|       1.30|       0.95|       1.04|       0.98|       0.88|       1.17|       0.83|       1.12|       1.14|       0.88|       0.89|       0.92|
|RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)                 |          1|       0.92|       0.92|       0.92|       0.99|       0.98|       0.86|         NA|       0.95|       0.82|       0.94|       0.81|       0.91|       0.91|       0.82|       0.91|       0.91|
|RollingWindow::RollingSum(x, wins, na_method = "ignore")                     |          1|       0.83|       0.93|       0.93|       1.04|       0.95|       0.96|         NA|       0.86|       0.97|       1.08|       0.94|       0.97|       0.93|       0.84|       0.87|       0.91|
|unjoin_cent_sums(rs3, rs2)                                                   |          1|       0.86|       0.75|       0.75|       0.87|       1.05|       0.95|       0.77|       0.78|       0.82|       0.84|       0.73|       0.81|       0.99|       0.78|       0.92|       0.91|
|roll::roll_mean(xm, wins)                                                    |          1|       0.95|       0.91|       0.91|       1.01|       1.00|       0.88|         NA|       0.91|       0.81|       0.96|       0.84|       0.93|       0.92|       0.83|       0.94|       0.90|
|as.centsums(x1, 4)                                                           |          1|       0.84|       0.86|       0.86|       0.92|       0.99|       0.91|       0.89|       0.93|       0.86|       0.90|       0.84|       0.95|       0.95|       0.96|       0.95|       0.88|

```r
allt %>% select(-sernum, -relchange, -meantime, -sumx_time, 
    -normalized) %>% distinct(expr, .keep_all = TRUE) %>% 
    arrange(desc(last_status)) %>% head(n = 20) %>% 
    kable()
```



|expr                                                                         | first_norm| last_norm| last_status|
|:----------------------------------------------------------------------------|----------:|---------:|-----------:|
|running_sum(x, wins)                                                         |       0.51|      1.40|         2.8|
|running_sd3(x, wins)                                                         |       8.72|     18.05|         2.1|
|sd3(x)                                                                       |      10.52|     14.00|         1.3|
|join_cent_sums(rs1, rs2)                                                     |       0.03|      0.04|         1.2|
|running_skew4(x, wins)                                                       |      48.07|     55.90|         1.2|
|running_tstat(x, wins)                                                       |       9.32|     10.82|         1.2|
|running_skew(x, wins)                                                        |      39.66|     45.79|         1.1|
|running_scaled(x, wins)                                                      |       9.03|     10.41|         1.1|
|running_zscored(x, wins)                                                     |       9.49|     10.65|         1.1|
|running_sharpe(x, wins)                                                      |       9.22|     10.34|         1.1|
|kurt5(x)                                                                     |     173.12|    193.52|         1.1|
|running_kurt5(x, wins)                                                       |      76.25|     84.85|         1.1|
|skew4(x)                                                                     |      97.66|    108.53|         1.1|
|unjoin_cent_sums(rs3, rs1)                                                   |       0.03|      0.03|         1.1|
|sd3(x, wts = w)                                                              |      13.06|     14.13|         1.1|
|kurtosis(x)                                                                  |      98.46|    102.81|         1.0|
|running_apx_median(x, wins)                                                  |     179.28|    184.18|         1.0|
|skewness(x)                                                                  |     100.16|    102.68|         1.0|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |     200.65|    201.66|         1.0|
|dumbk(x)                                                                     |     209.85|    210.87|         1.0|


