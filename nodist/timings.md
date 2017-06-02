

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
##                                "c5310a787373"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.4000      
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
##       sum(x)    80    80    82     80    81   111   100 a     
##      mean(x)   162   163   171    165   168   331   100 a     
##        sd(x)   459   463   484    474   490   602   100  b    
##  skewness(x)  8271  8367  8702   8464  8693 11799   100    d  
##  kurtosis(x)  8115  8201  8532   8271  8492 10293   100    d  
##       sd3(x)   855   859   892    872   902  1237   100   c   
##     skew4(x)  8334  8414  8608   8470  8673 10591   100    d  
##     kurt5(x) 14969 15102 15420  15246 15467 18538   100     e 
##     dumbk(x) 17036 17289 17879  17503 18642 20026   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16152 16342 16642  16472 16696 19933   100   c
##                                                               sd3(x, wts = w)   981   985  1026    999  1054  1318   100 a  
##                                                                 slow_sd(x, w)  1396  1415  1900   1511  2602  3740   100  b
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
##  as.centsums(x1, 1)  187  189  197    190  203  267   100  b  
##  as.centsums(x1, 2)  113  115  122    116  124  287   100 a   
##  as.centsums(x1, 3)  848  851  883    863  912 1069   100   c 
##  as.centsums(x1, 4) 1481 1484 1525   1507 1550 1683   100    d
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
##  c(obj1, obj2)  14 15   19     16 17 179   100   b
##  obj3 %-% obj1  11 12   14     12 13  94   100  a
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
##    join_cent_sums(rs1, rs2) 2.1 2.2  2.7    2.3 2.5  26   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.4    2.1 2.3  15   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.2    2.1 2.3   6   100   a
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
##  as.centcosums(x1, max_ord)  51 53   59     54 60 152   100   b
##             mobj3 %-% mobj1  16 18   21     19 20  45   100  a
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
    running_sd3(x, wins), running_skew4(x, wins), running_kurt5(x, 
        wins), running_tstat(x, wins), running_zscored(x, 
        wins), running_sharpe(x, wins), running_apx_median(x, 
        wins), running_centered(x, wins), running_scaled(x, 
        wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr   min    lq  mean median    uq    max neval        cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 31890 33351 37102  34772 36851 115571   100         i 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64592 67273 73110  69490 73426 163155   100          j
##                                           running_sum(x, wins)    72    77   113     83    95   2622   100 a         
##                                          running_mean(x, wins)    71    77    89     83    95    162   100 a         
##                                       roll::roll_sum(xm, wins)  1889  1951  2142   2021  2163   4899   100   c e     
##                                      roll::roll_mean(xm, wins)  2077  2132  2344   2232  2486   3897   100    de     
##                                        roll::roll_sd(xm, wins)  5728  5820  6216   5970  6491   8902   100      f    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   349   388   593    432   472   3617   100 a cd      
##                             RollingWindow::RollingSum(x, wins)   107   122   213    135   181   2053   100 a         
##                            RollingWindow::RollingMean(x, wins)   140   157   223    191   228   2130   100 ab        
##                             RollingWindow::RollingStd(x, wins)   223   253   367    304   336   2471   100 a c       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1743  1950  2118   2007  2168   4182   100  bcd      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1718  1976  2167   2064  2216   4662   100   c e     
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10165 10547 11533  10937 12340  16404   100       g   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   379   390   426    414   441    686   100 a c       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   440   452   516    473   517   2475   100 a cd      
##                                           running_sd3(x, wins)   612   632   716    653   684   3697   100 a cd      
##                                         running_skew4(x, wins)  3620  3697  4025   3777  4062   7023   100     e     
##                                         running_kurt5(x, wins)  5844  5959  6345   6057  6605   8276   100      f    
##                                         running_tstat(x, wins)   659   688   758    718   765   2538   100 a cd      
##                                       running_zscored(x, wins)   664   688   785    705   753   2792   100 a cd      
##                                        running_sharpe(x, wins)   663   678   744    711   768   1086   100 a cd      
##                                    running_apx_median(x, wins) 13884 14268 15594  14987 16493  23323   100        h  
##                                      running_centered(x, wins)   563   586   635    619   652   1115   100 a cd      
##                                        running_scaled(x, wins)   665   680   729    702   738   1245   100 a cd
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

