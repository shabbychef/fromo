

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
##                                "b0d7611d0fea"                                      "x86_64"                                     "unknown" 
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
##       sum(x)    79    80    88     81    91   174   100 a     
##      mean(x)   162   165   186    168   190   392   100 ab    
##        sd(x)   458   471   512    484   513   841   100  b    
##  skewness(x)  8305  8456  8944   8541  9075 12414   100    d  
##  kurtosis(x)  8157  8293  8754   8403  8975 10972   100    d  
##       sd3(x)   854   863   916    882   951  1166   100   c   
##     skew4(x)  8310  8477  8883   8546  9207 11323   100    d  
##     kurt5(x) 15121 15313 15974  15551 16136 19305   100     e 
##     dumbk(x) 17229 17403 18622  17980 19173 24284   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 16221 16602 17730  16973 18772 24793   100   c
##                                                               sd3(x, wts = w)   981   996  1058   1017  1088  1745   100 a  
##                                                                 slow_sd(x, w)  1424  1490  2225   1738  3056  4775   100  b
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
##  as.centsums(x1, 1)  187  188  197    190  195  374   100  b  
##  as.centsums(x1, 2)  113  114  121    115  124  288   100 a   
##  as.centsums(x1, 3)  847  850  873    863  876 1209   100   c 
##  as.centsums(x1, 4) 1470 1475 1509   1490 1533 1752   100    d
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
##  c(obj1, obj2)  14 15   20     16 16 232   100   a
##  obj3 %-% obj1  11 12   15     12 13 101   100   a
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
##    join_cent_sums(rs1, rs2) 2.1 2.3  3.3    2.4 2.5  32   100   a
##  unjoin_cent_sums(rs3, rs2) 1.9 2.0  2.4    2.2 2.4  11   100   a
##  unjoin_cent_sums(rs3, rs1) 1.9 2.0  2.8    2.1 2.4  41   100   a
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
##  as.centcosums(x1, max_ord)  53 55   63     56 69 202   100   b
##             mobj3 %-% mobj1  17 18   21     19 20  48   100  a
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
        n = wins, align = "right", fill = NA), running_sd3(x, 
        wins), running_skew4(x, wins), running_kurt5(x, 
        wins), running_tstat(x, wins), running_zscored(x, 
        wins), running_sharpe(x, wins), running_apx_median(x, 
        wins), running_centered(x, wins), running_scaled(x, 
        wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr   min    lq  mean median    uq    max neval      cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 32487 33848 37681  35457 40161  51696   100       g 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 65884 68277 77476  69902 80161 191345   100        h
##                                           running_sum(x, wins)    74    79   111     86   101   1955   100 a       
##                                          running_mean(x, wins)    71    78    90     88    98    126   100 a       
##                                       roll::roll_sum(xm, wins)  1800  1865  2064   1921  2211   4068   100 ab      
##                                      roll::roll_mean(xm, wins)  1975  2047  2245   2107  2317   3375   100 ab      
##                                        roll::roll_sd(xm, wins)  5443  5553  6102   5720  6592   8878   100   cd    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   361   396   571    439   520   2828   100 a       
##                             RollingWindow::RollingSum(x, wins)   113   138   221    159   213   1998   100 a       
##                            RollingWindow::RollingMean(x, wins)   143   170   246    209   246   2311   100 a       
##                             RollingWindow::RollingStd(x, wins)   233   263   377    301   339   2930   100 a       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1766  1938  2110   2007  2191   4385   100 ab      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1739  1986  2219   2049  2140   6672   100 ab      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10212 10788 12142  11474 13221  19346   100     e   
##                                           running_sd3(x, wins)   613   625   772    647   744   3431   100 a       
##                                         running_skew4(x, wins)  3639  3738  4056   3826  4326   6322   100  bc     
##                                         running_kurt5(x, wins)  5919  6054  6642   6158  7118  10401   100    d    
##                                         running_tstat(x, wins)   659   673   737    694   804    988   100 a       
##                                       running_zscored(x, wins)   663   682   795    714   797   3043   100 a       
##                                        running_sharpe(x, wins)   658   674   746    699   782   2346   100 a       
##                                    running_apx_median(x, wins) 13986 14184 15352  14392 16297  22449   100      f  
##                                      running_centered(x, wins)   573   581   673    603   677   1448   100 a       
##                                        running_scaled(x, wins)   658   678   800    725   826   2949   100 a
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

