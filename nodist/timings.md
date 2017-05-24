

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

print(sessionInfo())
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.2 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] methods stats   utils   base   
## 
## other attached packages:
##  [1] RcppParallel_4.3.20    RcppRoll_0.2.2         roll_1.0.7             RollingWindow_0.2      microbenchmark_1.4-2.1
##  [6] moments_0.14           dplyr_0.5.0            fromo_0.1.3.3001       ggplot2_2.2.1          knitr_1.15.1          
## [11] devtools_1.12.0        Quandl_2.8.0           xts_0.9-7              zoo_1.7-12             drat_0.1.2            
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.10     magrittr_1.5     grDevices_3.3.2  munsell_0.4.3    colorspace_1.3-2 lattice_0.20-33  R6_2.2.0        
##  [8] stringr_1.2.0    httr_1.2.1       plyr_1.8.4       tools_3.3.2      grid_3.3.2       gtable_0.2.0     DBI_0.6-1       
## [15] withr_1.0.2      assertthat_0.2.0 lazyeval_0.2.0   digest_0.6.12    tibble_1.3.0     formatR_1.5      graphics_3.3.2  
## [22] memoise_1.1.0    evaluate_0.10    stringi_1.1.5    scales_0.4.1     jsonlite_1.4
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
##       sum(x)    80    80    85     80    81   173   100 a    
##      mean(x)   162   164   175    166   176   278   100 a    
##        sd(x)   458   464   493    476   495   691   100 ab   
##  skewness(x)  7944  8067  8474   8171  8572 10652   100   c  
##  kurtosis(x)  7866  8037  8330   8103  8237 10074   100   c  
##       sd3(x)   854   859   890    870   909  1046   100  b   
##     skew4(x)  7957  8069  8262   8147  8299 10014   100   c  
##     kurt5(x) 14241 14406 14645  14487 14648 16957   100    d 
##     dumbk(x) 16549 16738 17753  16920 18085 56184   100     e
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
##                                                                          expr   min    lq  mean median    uq   max neval
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 15635 16107 16975  16622 17395 20612   100
##                                                               sd3(x, wts = w)   982  1026  1105   1079  1148  1519   100
##                                                                 slow_sd(x, w)  1391  1557  4171   2175  2906 95556   100
##  cld
##    c
##  a  
##   b
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
##  as.centsums(x1, 1)  198  203  240    219  250  523   100  b  
##  as.centsums(x1, 2)  121  130  155    140  156  457   100 a   
##  as.centsums(x1, 3)  831  899  998    953 1036 2014   100   c 
##  as.centsums(x1, 4) 1489 1563 1709   1647 1782 3266   100    d
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
##  c(obj1, obj2)  14 16   21     16 18 225   100   b
##  obj3 %-% obj1  11 13   16     13 14 102   100  a
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
##    join_cent_sums(rs1, rs2) 2.0 2.2  2.7    2.3 2.6  25   100   a
##  unjoin_cent_sums(rs3, rs2) 1.8 2.0  2.8    2.2 2.3  38   100   a
##  unjoin_cent_sums(rs3, rs1) 1.8 2.0  2.5    2.1 2.3  25   100   a
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
##  as.centcosums(x1, max_ord)  50 54   64     56 67 188   100   b
##             mobj3 %-% mobj1  16 18   22     19 21  57   100  a
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

checkit <- microbenchmark(sum(x), silly_fun(x, wins, 
    sum, na.rm = FALSE), silly_fun(x, wins, mean, na.rm = FALSE), 
    running_sum(x, wins), running_mean(x, wins), roll::roll_sum(xm, 
        wins), roll::roll_mean(xm, wins), roll::roll_sd(xm, 
        wins), RollingWindow::RollingSum(x, wins, na_method = "ignore"), 
    RollingWindow::RollingSum(x, wins), RollingWindow::RollingMean(x, 
        wins), RcppRoll::roll_sum(xm, n = wins, align = "right", 
        fill = NA), RcppRoll::roll_mean(xm, n = wins, 
        align = "right", fill = NA), running_skew4(x, 
        wins), running_kurt5(x, wins), running_tstat(x, 
        wins), running_zscored(x, wins), running_sharpe(x, 
        wins), running_apx_median(x, wins), running_centered(x, 
        wins), running_scaled(x, wins))
print(checkit)
```

```
## Unit: microseconds
##                                                           expr     min      lq  mean median    uq    max neval      cld
##                                                         sum(x)     8.3     9.7    12     11    13     36   100 a       
##                         silly_fun(x, wins, sum, na.rm = FALSE) 33860.9 36607.2 42388  38990 43894 123861   100       g 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 69576.1 73687.3 79847  77732 83472 148432   100        h
##                                           running_sum(x, wins)    23.7    30.8    43     36    49     96   100 a       
##                                          running_mean(x, wins)   146.2   159.2   200    168   197   2090   100 a       
##                                       roll::roll_sum(xm, wins)  1921.1  2020.9  2308   2088  2429   6508   100   cd    
##                                      roll::roll_mean(xm, wins)  2030.7  2165.4  2402   2287  2464   3596   100   cd    
##                                        roll::roll_sd(xm, wins)  5844.8  6007.8  6577   6223  6851  10490   100     e   
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   351.0   405.5   586    460   550   2852   100 a c     
##                             RollingWindow::RollingSum(x, wins)   109.6   138.6   314    193   226   2955   100 ab      
##                            RollingWindow::RollingMean(x, wins)   143.6   175.8   266    212   249   2920   100 a       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1723.4  2013.1  2208   2113  2366   4960   100  bcd    
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1752.4  2060.8  2268   2188  2373   5624   100   cd    
##                                         running_skew4(x, wins)  3451.0  3682.8  4066   3933  4314   6393   100    d    
##                                         running_kurt5(x, wins)  5594.8  5855.4  6450   6153  6880  10338   100     e   
##                                         running_tstat(x, wins)   665.0   707.9   788    770   846   1102   100 a c     
##                                       running_zscored(x, wins)   665.6   711.7   803    748   837   3151   100 a c     
##                                        running_sharpe(x, wins)   664.6   704.6   780    731   791   2989   100 a c     
##                                    running_apx_median(x, wins) 13195.0 13828.1 15166  14697 16108  24403   100      f  
##                                      running_centered(x, wins)   546.5   588.2   649    621   702    871   100 a c     
##                                        running_scaled(x, wins)   662.5   697.3   764    732   813   1041   100 a c
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

