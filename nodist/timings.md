

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
##                                        sysname                                        release                                        version 
##                                        "Linux"                            "4.4.0-137-generic" "#163-Ubuntu SMP Mon Sep 24 13:14:43 UTC 2018" 
##                                       nodename                                        machine                                          login 
##                                 "d5f3a1c82398"                                       "x86_64"                                      "unknown" 
##                                           user                                 effective_user 
##                                       "docker"                                       "docker"
```

```r
print(sessionInfo())
```

```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux buster/sid
## 
## Matrix products: default
## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
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
##  [1] RcppParallel_4.4.2   RcppRoll_0.3.0       roll_1.1.1           RollingWindow_0.2    microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8         
##  [8] ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6668    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       bindr_0.1.1      magrittr_1.5     grDevices_3.5.2  tidyselect_0.2.5 munsell_0.5.0    colorspace_1.3-2 R6_2.3.0         rlang_0.3.0.1   
## [10] stringr_1.3.1    plyr_1.8.4       tools_3.5.2      grid_3.5.2       gtable_0.2.0     xfun_0.4         withr_2.1.2      stats_3.5.2      lazyeval_0.2.1  
## [19] assertthat_0.2.0 tibble_1.4.2     crayon_1.3.4     bindrcpp_0.2.2   formatR_1.5      purrr_0.2.5      graphics_3.5.2   glue_1.3.0       evaluate_0.12   
## [28] stringi_1.2.4    compiler_3.5.2   pillar_1.3.1     scales_1.0.0     pkgconfig_2.0.2
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
##       sum(x)    76    77    84     79    85   136   100 a     
##      mean(x)   154   157   171    167   175   397   100 a     
##        sd(x)   312   345   464    359   396  9776   100 ab    
##  skewness(x)  6956  7191  7500   7333  7477 11276   100   c   
##  kurtosis(x)  6844  7147  7541   7279  7430 20609   100   c   
##       sd3(x)   660   673   703    692   729   815   100  b    
##     skew4(x)  7621  7923  8142   8029  8192 12127   100    d  
##     kurt5(x) 13408 13675 14043  13885 14052 20347   100     e 
##     dumbk(x) 14648 15014 15616  15181 15405 28698   100      f
```

```r
resdf <- checkit
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1746053   93    3.5e+06  186  2.3e+06  124
## Vcells 3586525   27    2.2e+07  168  3.4e+07  262
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 13783 13990 14182  14141 14298 15055   100   c
##                                                               sd3(x, wts = w)   798   809   852    852   882  1004   100 a  
##                                                                 slow_sd(x, w)   867   976  1654   1694  1785  6527   100  b
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1747469   93    3.5e+06  186  2.3e+06  124
## Vcells 3592615   28    1.4e+07  107  3.4e+07  262
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
##  as.centsums(x1, 1)   73   75   79     76   83  103   100 a   
##  as.centsums(x1, 2)  109  111  123    114  124  386   100  b  
##  as.centsums(x1, 3)  762  769  819    812  854  946   100   c 
##  as.centsums(x1, 4) 1324 1353 1424   1424 1477 1621   100    d
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1747954   93    3.5e+06  186  2.3e+06  124
## Vcells 3616379   28    1.4e+07  107  3.4e+07  262
```

```r
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
##  c(obj1, obj2)  14 15   19     16 17 216   100   a
##  obj3 %-% obj1  11 12   14     13 13 112   100   a
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1748633   93    3.5e+06  186  2.3e+06  124
## Vcells 3617668   28    1.4e+07  107  3.4e+07  262
```

```r
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
##    join_cent_sums(rs1, rs2) 2.1 2.2  2.8    2.3 2.6 34.7   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.1  2.3    2.2 2.5  3.3   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.1  2.5    2.3 2.6 10.7   100   a
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1748712   93    3.5e+06  186  2.3e+06  124
## Vcells 3618431   28    1.4e+07  107  3.4e+07  262
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
##  as.centcosums(x1, max_ord)  52 54   58     55 56 209   100   b
##             mobj3 %-% mobj1  16 17   19     18 19  45   100  a
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1748831   93    3.5e+06  186  2.3e+06  124
## Vcells 3599093   28    1.4e+07  107  3.4e+07  262
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 17462 18036 22117  18383 24733 125325   100        h 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 45834 47521 50510  48584 53814  60639   100         i
##                                           running_sum(x, wins)   100   107   120    118   129    160   100 a        
##                                          running_mean(x, wins)   103   107   117    112   121    199   100 a        
##                                       roll::roll_sum(xm, wins)    90   116   141    138   163    326   100 a        
##                                      roll::roll_mean(xm, wins)    97   121   155    145   181    597   100 a        
##                                        roll::roll_sd(xm, wins)   368   413   454    444   489    610   100 a        
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   240   278   382    311   345   7197   100 a        
##                             RollingWindow::RollingSum(x, wins)   104   134   171    162   206    308   100 a        
##                            RollingWindow::RollingMean(x, wins)   133   169   261    206   225   5968   100 a        
##                             RollingWindow::RollingStd(x, wins)   217   248   363    292   322   7584   100 a        
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1630  1779  1920   1860  1931   8527   100   c      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1666  1779  1880   1858  1964   2510   100  bc      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA)  8062  8508  9392   8735  8996  15642   100      f   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   559   601   661    632   681   1223   100 a        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   588   629   680    668   704   1153   100 ab       
##                                           running_sd3(x, wins)   745   777   898    824   863   6098   100 a c      
##                                          running_skew(x, wins)  3601  3750  3873   3820  3976   4821   100    d     
##                                         running_skew4(x, wins)  3647  3805  3927   3913  4016   4323   100    d     
##                                          running_kurt(x, wins)  5134  5308  5458   5430  5570   6120   100     e    
##                                         running_kurt5(x, wins)  5652  5881  6001   5984  6155   6373   100     e    
##                                         running_tstat(x, wins)   779   820   940    864   898   6421   100 a c      
##                                       running_zscored(x, wins)   742   788   893    826   875   6115   100 a c      
##                                        running_sharpe(x, wins)   749   782   916    828   895   7364   100 a c      
##                                    running_apx_median(x, wins) 15314 15581 15915  15777 15970  21530   100       g  
##                                      running_centered(x, wins)  1001  1086  1158   1145  1208   2139   100 a c      
##                                        running_scaled(x, wins)   759   804   844    829   878   1090   100 a c
```

```r
resdf <- rbind(resdf, checkit)
gc()
```

```
##           used (Mb) gc trigger (Mb) max used (Mb)
## Ncells 1756743   94    3.5e+06  186  3.5e+06  186
## Vcells 3179434   24    1.4e+07  107  3.4e+07  262
```


```r
library(readr)
FAKE_IT <- FALSE
if (FAKE_IT) {
    resdf <- readr::read_csv("timings.csv")
    # print(resdf)
} else {
    readr::write_csv(resdf, "timings.csv")
}
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
    mutate(relchange = normalized/first_norm, last_status = last_norm/first_norm) %>% 
    mutate(exgrp = ifelse(grepl("^(roll|RollingWindow|RcppRoll)::", 
        expr), "brand_x", ifelse(grepl("^running_", 
        expr), "running", "summarizing")))

library(ggplot2)
ph <- allt %>% dplyr::filter(!grepl("brand_x", exgrp)) %>% 
    ggplot(aes(sernum, normalized, group = expr, color = expr)) + 
    geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + facet_grid(exgrp ~ ., scales = "free") + 
    labs(x = "release", y = "mean time taken, relative to sum(x)", 
        title = "fromo microbenchmark timings, lasagna")
print(ph)
```

<img src="figure/all_timing_stats-1.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="700px" height="600px" />

```r
ph <- allt %>% dplyr::filter(!grepl("brand_x", exgrp)) %>% 
    ggplot(aes(sernum, relchange, group = expr, color = expr)) + 
    geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + facet_grid(exgrp ~ ., scales = "free") + 
    labs(x = "release", y = "normalized time taken, relative to first iteration", 
        title = "fromo microbenchmark timings, spaghetti")
print(ph)
```

<img src="figure/all_timing_stats-2.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="700px" height="600px" />

```r
ph <- allt %>% dplyr::filter(!grepl("brand_x", exgrp)) %>% 
    ggplot(aes(sernum, relchange)) + geom_boxplot(aes(group = sernum), 
    alpha = 0.7) + stat_summary(aes(group = "1", color = "mean"), 
    fun.y = mean, geom = "line") + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    facet_grid(exgrp ~ ., scales = "free") + labs(x = "release", 
    y = "normalized time taken, relative to first iteration", 
    color = "stat", title = "fromo microbenchmark timings, boxplots")
print(ph)
```

<img src="figure/all_timing_stats-3.png" title="plot of chunk all_timing_stats" alt="plot of chunk all_timing_stats" width="700px" height="600px" />

```r
allt %>% dplyr::filter(!grepl("brand_x", exgrp)) %>% 
    select(expr, sernum, relchange, last_status) %>% 
    tidyr::spread(key = "sernum", value = "relchange") %>% 
    arrange(desc(last_status)) %>% select(-last_status) %>% 
    head(n = 50) %>% kable()
```



|expr                                                                         | 0.1.3.3001| 0.1.3.3200| 0.1.3.3330| 0.1.3.3400| 0.1.3.4000| 0.1.3.4100| 0.1.3.4200| 0.1.3.4300| 0.1.3.4400| 0.1.3.4500| 0.1.3.4510| 0.1.3.5000| 0.1.3.5010| 0.1.3.5020| 0.1.3.6000| 0.1.3.6100| 0.1.3.6200| 0.1.3.6300| 0.1.3.6500| 0.1.3.6510| 0.1.3.6520| 0.1.3.6530| 0.1.3.6540| 0.1.3.6550| 0.1.3.6560| 0.1.3.6600| 0.1.3.6660| 0.1.3.6661| 0.1.3.6662| 0.1.3.6665| 0.1.3.6666| 0.1.3.6667| 0.1.3.6668|
|:----------------------------------------------------------------------------|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
|running_sum(x, wins)                                                         |          1|         NA|         NA|         NA|       2.72|       2.18|       1.95|         NA|       2.02|       2.53|       2.74|       2.45|       3.23|       3.72|       2.50|       2.70|       2.74|       2.74|       2.79|       2.74|       2.76|       2.76|       2.84|       2.70|       2.77|       3.31|       2.76|       2.74|       3.29|       3.14|       2.80|       3.10|       3.00|
|running_centered(x, wins)                                                    |          1|       0.93|       0.91|       0.91|       1.01|       0.98|       0.92|         NA|       0.94|       0.84|       0.91|       0.87|       1.08|       1.06|       0.94|       1.04|       0.95|       0.95|       1.56|       1.61|       1.94|       1.86|       2.03|       2.05|       2.02|       2.04|       1.96|       1.89|       1.82|       1.94|       2.03|       1.83|       1.78|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)                  |         NA|         NA|         NA|         NA|       1.00|       1.03|       0.86|         NA|       0.90|       0.85|       0.92|       0.83|       1.07|       1.06|       1.07|       1.02|       1.00|       1.00|       1.30|       1.35|       1.94|       2.04|       1.90|       1.90|       1.77|       1.74|       1.73|       1.81|       1.50|       1.71|       1.43|       1.55|       1.54|
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)                    |         NA|         NA|         NA|         NA|       1.00|       1.03|       0.87|         NA|       0.87|       0.76|       0.87|       0.75|       1.17|       1.08|       1.07|       0.96|       0.89|       0.89|       1.20|       1.24|       1.69|       1.79|       1.71|       1.64|       1.51|       1.52|       1.48|       1.61|       1.24|       1.41|       1.20|       1.31|       1.32|
|running_sd3(x, wins)                                                         |         NA|         NA|         NA|         NA|       1.00|       1.04|       0.86|         NA|       0.95|       0.78|       0.88|       0.80|       2.14|       2.03|       1.88|       2.05|       2.03|       2.03|       1.56|       1.46|       1.85|       1.60|       1.76|       1.69|       1.44|       1.59|       1.43|       1.35|       1.47|       1.47|       1.22|       1.33|       1.26|
|running_tstat(x, wins)                                                       |          1|       0.92|       0.94|       0.94|       0.99|       1.01|       0.82|         NA|       0.86|       0.80|       0.90|       0.76|       1.23|       1.23|       1.10|       1.14|       1.18|       1.18|       1.13|       1.13|       1.43|       1.50|       1.48|       1.44|       1.33|       1.28|       1.28|       1.35|       1.19|       1.38|       1.16|       1.31|       1.22|
|running_sharpe(x, wins)                                                      |          1|       0.99|       0.90|       0.90|       0.98|       1.02|       0.84|         NA|       0.90|       0.79|       0.88|       0.77|       1.20|       1.19|       1.05|       1.11|       1.13|       1.13|       1.11|       1.08|       1.55|       1.42|       1.40|       1.44|       1.28|       1.31|       1.29|       1.34|       1.12|       1.29|       1.11|       1.26|       1.19|
|running_skew(x, wins)                                                        |         NA|         NA|         NA|         NA|         NA|       1.00|       1.09|         NA|       1.16|       1.05|       1.18|       1.04|       1.17|       1.18|       1.05|       1.18|       1.13|       1.13|       1.29|       1.29|       1.36|       1.35|       1.39|       1.38|       1.34|       1.39|       1.93|       1.93|       1.72|       1.75|       1.70|       1.19|       1.17|
|running_scaled(x, wins)                                                      |          1|       0.94|       0.91|       0.91|       0.98|       1.03|       0.88|         NA|       0.92|       0.79|       0.88|       0.79|       1.19|       1.19|       1.08|       1.19|       1.15|       1.15|       1.15|       1.15|       1.53|       1.47|       1.48|       1.53|       1.27|       1.40|       1.31|       1.31|       1.14|       1.25|       1.09|       1.20|       1.14|
|running_zscored(x, wins)                                                     |          1|       0.90|       0.88|       0.88|       1.01|       0.98|       0.80|         NA|       0.89|       0.77|       0.84|       0.78|       1.20|       1.17|       1.02|       1.14|       1.10|       1.10|       1.15|       1.09|       1.47|       1.41|       1.44|       1.42|       1.28|       1.29|       1.25|       1.26|       1.09|       1.29|       1.07|       1.20|       1.13|
|running_apx_median(x, wins)                                                  |          1|       0.92|       0.90|       0.90|       1.06|       1.09|       0.93|         NA|       1.00|       0.89|       1.00|       0.88|       1.05|       1.05|       0.93|       1.05|       1.02|       1.02|       1.04|       1.00|       1.02|       1.03|       1.01|       1.01|       1.01|       1.04|       1.01|       1.01|       1.01|       1.08|       0.99|       1.08|       1.07|
|join_cent_sums(rs1, rs2)                                                     |          1|       1.12|       0.96|       0.96|       1.02|       1.25|       1.12|       0.97|       1.10|       1.03|       1.05|       0.93|       1.02|       1.14|       1.01|       1.11|       1.16|       1.16|       1.11|       1.03|       1.13|       1.04|       1.06|       1.25|       1.05|       1.08|       1.15|       1.05|       1.16|       0.97|       1.18|       1.02|       1.04|
|skew4(x)                                                                     |          1|       1.01|       1.01|       1.01|       1.07|       1.07|       1.06|       1.05|       1.10|       1.08|       1.07|       1.09|       1.14|       1.14|       1.12|       1.13|       1.09|       1.09|       1.16|       1.14|       1.13|       1.13|       1.15|       1.16|       1.12|       1.15|       1.77|       1.76|       1.59|       1.55|       1.51|       1.01|       1.00|
|sum(x)                                                                       |          1|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|       1.00|
|mean(x)                                                                      |          1|       1.02|       1.00|       1.00|       1.01|       1.01|       1.01|       1.01|       1.01|       1.06|       0.97|       1.03|       1.02|       1.00|       1.05|       1.02|       1.00|       1.00|       1.00|       0.98|       0.98|       0.99|       1.02|       1.01|       0.98|       1.01|       0.97|       1.02|       1.00|       1.01|       1.03|       0.97|       0.99|
|running_skew4(x, wins)                                                       |          1|       0.94|       0.89|       0.89|       1.02|       1.03|       0.92|         NA|       0.97|       0.87|       0.96|       0.86|       1.17|       1.18|       1.04|       1.17|       1.13|       1.13|       1.15|       1.13|       1.16|       1.19|       1.19|       1.16|       1.13|       1.17|       1.62|       1.61|       1.41|       1.45|       1.40|       0.99|       0.99|
|unjoin_cent_sums(rs3, rs1)                                                   |          1|       1.21|       0.97|       0.97|       0.91|       1.34|       0.88|       0.99|       0.97|       0.99|       0.98|       0.85|       1.10|       1.19|       0.99|       1.11|       1.09|       1.09|       1.08|       1.04|       0.94|       1.05|       0.92|       1.20|       1.12|       1.00|       1.03|       0.90|       1.08|       0.99|       1.22|       1.00|       0.98|
|kurt5(x)                                                                     |          1|       1.02|       1.00|       1.00|       1.09|       1.10|       1.07|       1.06|       1.12|       1.09|       1.08|       1.08|       1.15|       1.16|       1.12|       1.14|       1.09|       1.09|       1.15|       1.12|       1.11|       1.12|       1.11|       1.14|       1.12|       1.14|       1.05|       1.04|       0.95|       0.92|       0.88|       0.98|       0.97|
|running_kurt5(x, wins)                                                       |          1|       0.94|       0.91|       0.91|       1.01|       1.03|       0.93|         NA|       1.00|       0.90|       0.98|       0.89|       1.13|       1.13|       1.01|       1.13|       1.09|       1.09|       1.13|       1.11|       1.11|       1.12|       1.09|       1.09|       1.10|       1.12|       1.05|       1.05|       0.96|       0.94|       0.92|       0.97|       0.95|
|running_kurt(x, wins)                                                        |         NA|         NA|         NA|         NA|         NA|       1.00|       0.91|         NA|       0.97|       0.85|       0.95|       0.86|       0.97|       0.96|       0.86|       0.96|       0.92|       0.92|       1.05|       1.02|       1.06|       1.08|       1.09|       1.08|       1.06|       1.09|       1.04|       1.03|       0.97|       0.95|       0.93|       0.95|       0.94|
|as.centcosums(x1, max_ord)                                                   |          1|       0.93|       1.02|       1.02|       0.96|       1.24|       0.98|       1.02|       0.98|       0.87|       1.13|       0.85|       1.12|       1.16|       0.83|       0.93|       0.92|       0.92|       1.06|       0.91|       1.07|       0.93|       0.94|       0.96|       0.98|       0.95|       1.00|       0.86|       0.87|       0.91|       0.91|       0.92|       0.92|
|kurtosis(x)                                                                  |          1|       1.00|       1.00|       1.00|       1.06|       1.06|       1.03|       1.00|       1.05|       1.04|       1.03|       1.03|       1.04|       1.07|       1.04|       1.04|       1.04|       1.04|       1.05|       1.03|       1.01|       1.03|       1.06|       1.09|       1.05|       1.06|       1.05|       1.02|       0.92|       0.91|       0.90|       0.92|       0.91|
|skewness(x)                                                                  |          1|       1.00|       1.00|       1.00|       1.06|       1.07|       1.04|       1.01|       1.06|       1.05|       1.03|       1.03|       1.06|       1.08|       1.06|       1.05|       0.98|       0.98|       1.08|       1.03|       1.03|       1.04|       1.04|       1.06|       1.03|       1.06|       1.04|       1.03|       0.90|       0.89|       0.88|       0.90|       0.89|
|dumbk(x)                                                                     |          1|       1.03|       1.01|       1.01|       1.04|       1.04|       1.01|       0.99|       1.04|       1.02|       1.02|       1.02|       1.04|       1.04|       1.03|       1.03|       0.96|       0.96|       1.04|       1.02|       1.02|       1.02|       1.02|       1.04|       1.02|       1.05|       1.03|       1.02|       0.89|       0.89|       0.89|       0.89|       0.89|
|mobj3 %-% mobj1                                                              |          1|       0.93|       1.03|       1.03|       0.97|       1.30|       0.95|       1.04|       0.98|       0.88|       1.17|       0.83|       1.12|       1.14|       0.88|       0.89|       0.93|       0.93|       1.05|       0.91|       1.05|       0.92|       0.94|       0.99|       1.03|       0.93|       0.98|       0.84|       0.82|       0.89|       0.86|       0.86|       0.88|
|obj3 %-% obj1                                                                |          1|       1.53|       0.96|       0.96|       0.91|       0.97|       1.00|       0.87|       1.00|       0.85|       0.85|       0.79|       0.92|       0.95|       0.90|       0.90|       0.96|       0.96|       0.97|       1.02|       0.92|       0.91|       0.83|       0.90|       0.91|       0.98|       1.04|       0.86|       0.87|       0.89|       1.41|       0.88|       0.88|
|c(obj1, obj2)                                                                |          1|       1.29|       0.97|       0.97|       0.93|       0.98|       1.04|       0.87|       0.95|       0.85|       0.85|       0.81|       0.89|       0.88|       0.83|       0.96|       0.93|       0.93|       0.88|       1.01|       0.87|       0.91|       0.92|       0.90|       0.91|       1.01|       1.03|       0.94|       0.87|       0.89|       1.37|       0.90|       0.88|
|cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) |          1|       0.93|       0.94|       0.94|       1.01|       1.07|       0.96|       0.95|       1.02|       0.96|       1.06|       0.91|       1.04|       1.01|       0.97|       1.00|       0.96|       0.96|       1.02|       1.00|       1.01|       1.03|       0.98|       1.04|       1.01|       1.03|       1.03|       1.00|       0.90|       0.86|       0.79|       0.86|       0.86|
|sd(x)                                                                        |          1|       1.02|       1.01|       1.01|       1.01|       1.01|       1.03|       0.99|       1.01|       1.02|       1.00|       1.02|       1.02|       1.02|       1.04|       1.02|       0.69|       0.69|       1.01|       0.99|       0.99|       1.00|       1.30|       1.49|       0.99|       1.04|       1.00|       1.03|       0.77|       0.75|       0.78|       0.73|       0.85|
|as.centsums(x1, 4)                                                           |          1|       0.84|       0.86|       0.86|       0.92|       0.99|       0.91|       0.89|       0.93|       0.86|       0.90|       0.84|       0.95|       0.95|       0.96|       0.95|       0.92|       0.92|       0.96|       0.96|       0.97|       0.98|       1.00|       0.91|       0.97|       0.94|       0.93|       0.90|       0.82|       0.80|       0.70|       0.84|       0.85|
|as.centsums(x1, 3)                                                           |          1|       0.85|       0.86|       0.86|       0.91|       0.96|       0.89|       0.87|       0.92|       0.84|       0.89|       0.83|       0.95|       0.95|       0.96|       0.94|       0.91|       0.91|       0.99|       0.97|       0.99|       1.00|       1.00|       0.93|       0.97|       0.96|       1.54|       1.46|       1.35|       1.31|       1.17|       0.84|       0.83|
|unjoin_cent_sums(rs3, rs2)                                                   |          1|       0.86|       0.75|       0.75|       0.87|       1.05|       0.95|       0.77|       0.78|       0.82|       0.84|       0.73|       0.81|       0.99|       0.78|       0.92|       0.93|       0.93|       0.86|       0.86|       0.81|       0.87|       0.84|       1.03|       0.83|       0.88|       1.05|       0.82|       0.92|       0.84|       1.08|       0.87|       0.83|
|sd3(x)                                                                       |          1|       1.03|       1.01|       1.01|       1.03|       1.03|       1.02|       1.00|       1.02|       1.01|       0.99|       1.00|       1.77|       1.80|       1.71|       1.76|       1.34|       1.34|       1.33|       1.30|       1.31|       1.31|       1.31|       1.35|       1.01|       1.01|       1.00|       0.99|       1.00|       0.98|       0.97|       1.00|       0.81|
|as.centsums(x1, 2)                                                           |          1|       0.82|       0.81|       0.81|       0.81|       0.88|       0.80|       0.78|       0.81|       0.77|       0.78|       0.72|       1.38|       1.38|       1.41|       1.38|       0.98|       0.98|       0.99|       0.98|       0.99|       1.02|       0.99|       0.94|       0.83|       0.82|       0.83|       0.78|       0.79|       0.80|       0.71|       0.79|       0.79|
|sd3(x, wts = w)                                                              |          1|       0.94|       0.94|       0.94|       0.96|       0.97|       0.90|       0.90|       0.95|       0.90|       0.97|       0.83|       1.63|       1.56|       1.50|       1.54|       1.07|       1.07|       1.07|       1.05|       1.04|       1.06|       1.03|       1.10|       1.05|       0.97|       0.97|       0.92|       0.96|       0.94|       0.86|       0.93|       0.79|
|running_sum(x, wins, robust = FALSE)                                         |         NA|       1.00|       0.66|       0.66|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|
|silly_fun(x, wins, mean, na.rm = FALSE)                                      |          1|       0.91|       0.89|       0.89|       0.94|       0.94|       0.88|         NA|       0.92|       0.82|       0.90|       0.80|       0.92|       0.88|       0.79|       0.88|       0.76|       0.76|       0.88|       0.88|       0.87|       0.89|       0.88|       0.87|       0.89|       0.92|       0.86|       0.90|       0.66|       0.71|       0.59|       0.67|       0.65|
|running_mean(x, wins)                                                        |          1|       0.92|       0.59|       0.59|       0.46|       0.47|       0.53|         NA|       0.56|       0.66|       0.72|       0.52|       0.59|       0.59|       0.53|       0.70|       0.59|       0.59|       0.81|       0.77|       0.60|       0.71|       0.63|       0.70|       0.59|       0.62|       0.71|       0.60|       0.60|       0.64|       0.54|       0.60|       0.60|
|silly_fun(x, wins, sum, na.rm = FALSE)                                       |          1|       0.92|       0.90|       0.90|       0.90|       0.90|       0.82|         NA|       0.87|       0.77|       0.87|       0.76|       0.86|       0.88|       0.76|       0.83|       0.67|       0.67|       0.86|       0.84|       0.84|       0.85|       0.87|       0.86|       0.84|       0.86|       0.82|       0.86|       0.49|       0.60|       0.49|       0.54|       0.53|
|running_sum(x, wins, robust = TRUE)                                          |         NA|       1.00|       0.48|       0.48|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|         NA|
|slow_sd(x, w)                                                                |          1|       0.94|       0.89|       0.89|       0.47|       0.49|       0.48|       0.47|       0.48|       0.47|       0.55|       0.44|       0.48|       0.47|       0.49|       0.46|       0.54|       0.54|       0.48|       0.47|       0.46|       0.46|       0.45|       0.49|       0.46|       0.49|       0.47|       0.46|       0.42|       0.41|       0.39|       0.39|       0.39|
|as.centsums(x1, 1)                                                           |          1|       0.81|       0.82|       0.82|       0.85|       0.92|       0.85|       0.81|       0.87|       0.83|       0.85|       0.77|       0.87|       0.85|       0.93|       0.87|       0.35|       0.35|       0.36|       0.35|       0.36|       0.36|       0.36|       0.34|       0.35|       0.36|       0.35|       0.34|       0.34|       0.34|       0.31|       0.34|       0.33|

```r
allt %>% dplyr::filter(!grepl("brand_x", exgrp)) %>% 
    select(-sernum, -relchange, -meantime, -sumx_time, 
        -normalized) %>% distinct(expr, .keep_all = TRUE) %>% 
    arrange(desc(last_status)) %>% head(n = 20) %>% 
    kable()
```



|expr                                                        | first_norm| last_norm| last_status|exgrp       |
|:-----------------------------------------------------------|----------:|---------:|-----------:|:-----------|
|running_sum(x, wins)                                        |       0.51|      1.52|        3.00|running     |
|running_centered(x, wins)                                   |       7.68|     13.66|        1.78|running     |
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) |       5.19|      7.99|        1.54|running     |
|running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   |       6.29|      8.33|        1.32|running     |
|running_sd3(x, wins)                                        |       8.72|     10.95|        1.26|running     |
|running_tstat(x, wins)                                      |       9.32|     11.33|        1.22|running     |
|running_sharpe(x, wins)                                     |       9.22|     10.99|        1.19|running     |
|running_skew(x, wins)                                       |      39.66|     46.57|        1.17|running     |
|running_scaled(x, wins)                                     |       9.03|     10.27|        1.14|running     |
|running_zscored(x, wins)                                    |       9.49|     10.75|        1.13|running     |
|running_apx_median(x, wins)                                 |     179.28|    192.13|        1.07|running     |
|join_cent_sums(rs1, rs2)                                    |       0.03|      0.03|        1.04|summarizing |
|skew4(x)                                                    |      97.66|     97.67|        1.00|summarizing |
|sum(x)                                                      |       1.00|      1.00|        1.00|summarizing |
|mean(x)                                                     |       2.07|      2.05|        0.99|summarizing |
|running_skew4(x, wins)                                      |      48.07|     47.58|        0.99|running     |
|unjoin_cent_sums(rs3, rs1)                                  |       0.03|      0.03|        0.98|summarizing |
|kurt5(x)                                                    |     173.12|    168.48|        0.97|summarizing |
|running_kurt5(x, wins)                                      |      76.25|     72.69|        0.95|running     |
|running_kurt(x, wins)                                       |      70.61|     66.45|        0.94|running     |


