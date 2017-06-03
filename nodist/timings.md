

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
##                                "12287af8901f"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4100      
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
##       sum(x)    80    80    83     80    84   123   100 a     
##      mean(x)   162   164   174    166   175   268   100 a     
##        sd(x)   458   463   490    476   496   720   100  b    
##  skewness(x)  8270  8347  8908   8500  9306 11811   100    d  
##  kurtosis(x)  8121  8217  8675   8333  8816 12821   100    d  
##       sd3(x)   851   861   897    875   913  1163   100   c   
##     skew4(x)  8310  8415  8664   8486  8752 11555   100    d  
##     kurt5(x) 14972 15185 15775  15318 15838 20302   100     e 
##     dumbk(x) 17012 17284 18173  17600 18709 23022   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16176 16389 17857  16628 19654 22352   100   c
##                                                               sd3(x, wts = w)   981   985  1053   1000  1066  1439   100 a  
##                                                                 slow_sd(x, w)  1396  1415  1994   1643  2602  3964   100  b
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
##  as.centsums(x1, 1)  188  201  217    204  217  379   100  b  
##  as.centsums(x1, 2)  114  122  134    125  137  290   100 a   
##  as.centsums(x1, 3)  850  917  943    924  951 1309   100   c 
##  as.centsums(x1, 4) 1473 1600 1656   1628 1672 2040   100    d
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
##  c(obj1, obj2)  16 16   20     17 18 195   100   b
##  obj3 %-% obj1  12 13   15     13 14 100   100  a
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
##    join_cent_sums(rs1, rs2) 2.3 2.6  3.3    2.9 3.1  37   100   a
##  unjoin_cent_sums(rs3, rs2) 2.0 2.3  2.9    2.5 2.7  27   100   a
##  unjoin_cent_sums(rs3, rs1) 2.0 2.4  3.3    2.6 2.8  63   100   a
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
##  as.centcosums(x1, max_ord)  56 61   78     71 85 179   100   b
##             mobj3 %-% mobj1  18 21   28     23 32  99   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 31688 33286 37348  35079 38068 116862   100           k 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64495 67243 73763  68962 72252 159787   100            l
##                                           running_sum(x, wins)    73    78    92     83    99    179   100 a           
##                                          running_mean(x, wins)    74    79    92     84    93    203   100 a           
##                                       roll::roll_sum(xm, wins)  1878  1930  2120   1982  2128   5048   100  b de       
##                                      roll::roll_mean(xm, wins)  2069  2117  2354   2192  2469   3696   100    d f      
##                                        roll::roll_sd(xm, wins)  5722  5846  6624   6160  7014  10986   100        h    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   365   395   545    435   516   3425   100 a  d        
##                             RollingWindow::RollingSum(x, wins)   108   122   248    148   197   3152   100 ab          
##                            RollingWindow::RollingMean(x, wins)   140   157   269    193   237   2306   100 ab          
##                             RollingWindow::RollingStd(x, wins)   226   244   373    288   329   2994   100 abc         
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1751  1948  2116   2023  2214   2970   100  b de       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1722  2000  2260   2063  2318   5538   100   cd f      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10228 10726 12965  11575 13115  83536   100         i   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   389   400   445    424   452   1040   100 a  d        
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   448   462   541    493   535   2708   100 a  d        
##                                           running_sd3(x, wins)   612   624   757    650   701   2859   100 a  d        
##                                          running_skew(x, wins)  2978  3053  3297   3178  3367   4577   100     ef      
##                                         running_skew4(x, wins)  3622  3674  4107   3812  4185   8628   100      fg     
##                                          running_kurt(x, wins)  5195  5309  5870   5607  6087  10587   100       gh    
##                                         running_kurt5(x, wins)  5845  5939  6526   6182  6833  10290   100        h    
##                                         running_tstat(x, wins)   670   682   780    715   771   2437   100 a  d        
##                                       running_zscored(x, wins)   669   693   774    730   805   1368   100 a  d        
##                                        running_sharpe(x, wins)   671   680   780    719   781   2734   100 a  d        
##                                    running_apx_median(x, wins) 13950 14249 16176  15068 17838  24733   100          j  
##                                      running_centered(x, wins)   564   574   623    599   629    997   100 a  d        
##                                        running_scaled(x, wins)   668   678   771    704   771   2528   100 a  d
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
```

```
## Error in library(tidyr): there is no package called 'tidyr'
```

```r
library(knitr)

allt <- data.frame(fname = dir(".", "*.csv"), stringsAsFactors = FALSE) %>% 
    filter(grepl("^timings_\\d.+\\d+.csv$", fname)) %>% 
    group_by(fname) %>% mutate(tims = list(readr::read_csv(fname))) %>% 
    ungroup() %>% tidyr::unnest() %>% mutate(sernum = gsub("^timings_(.+).csv$", 
    "\\1", fname)) %>% select(-fname) %>% group_by(sernum, 
    expr) %>% summarize(meantime = mean(time, na.rm = TRUE)) %>% 
    ungroup() %>% mutate(is_numeraire = grepl("^sum\\(x\\)$", 
    expr)) %>% arrange(!is_numeraire) %>% group_by(sernum) %>% 
    mutate(numv = first(meantime)) %>% ungroup() %>% 
    mutate(normalized = meantime/numv) %>% arrange(sernum) %>% 
    group_by(expr) %>% mutate(first_norm = first(normalized)) %>% 
    ungroup() %>% mutate(relchange = normalized/first_norm)
```

```
## Error in grepl("^timings_\\d.+\\d+.csv$", fname): object 'fname' not found
```

```r
library(ggplot2)
ph <- allt %>% ggplot(aes(sernum, normalized, group = expr, 
    color = expr)) + geom_line() + geom_point() + scale_y_log10() + 
    labs(x = "release", y = "mean time taken, relative to sum(x)", 
        title = "fromo microbenchmark timings")
```

```
## Error in eval(expr, envir, enclos): object 'allt' not found
```

```r
print(ph)
```

```
## Error in print(ph): object 'ph' not found
```

```r
ph <- allt %>% ggplot(aes(sernum, relchange, group = expr, 
    color = expr)) + geom_line() + geom_point() + scale_y_log10() + 
    labs(x = "release", y = "normalized time taken, relative to first iteration", 
        title = "fromo microbenchmark timings")
```

```
## Error in eval(expr, envir, enclos): object 'allt' not found
```

```r
print(ph)
```

```
## Error in print(ph): object 'ph' not found
```

```r
allt %>% arrange(sernum) %>% group_by(expr) %>% mutate(perfo = mean(relchange, 
    na.rm = TRUE), first_mean = first(meantime), last_mean = last(meantime)) %>% 
    ungroup() %>% distinct(expr, .keep_all = TRUE) %>% 
    select(expr, first_mean, last_mean, perfo) %>% 
    mutate(rel_mean = last_mean/first_mean) %>% arrange(desc(perfo)) %>% 
    head(n = 10) %>% kable()
```

```
## Error in eval(expr, envir, enclos): object 'allt' not found
```


