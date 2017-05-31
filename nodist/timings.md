

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
##  [6] moments_0.14           dplyr_0.5.0            fromo_0.1.3.3330       ggplot2_2.2.1          knitr_1.15.1          
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
##       sum(x)    80    80    84     80    84   126   100 a    
##      mean(x)   162   165   173    166   171   297   100 a    
##        sd(x)   458   467   494    481   506   634   100 ab   
##  skewness(x)  7922  8018  8366   8088  8236 10970   100   c  
##  kurtosis(x)  7779  7878  8271   7939  8265 10215   100   c  
##       sd3(x)   855   860   886    869   893  1056   100  b   
##     skew4(x)  7926  8046  8230   8101  8256  9777   100   c  
##     kurt5(x) 14088 14243 14449  14308 14404 19178   100    d 
##     dumbk(x) 16320 16555 17767  16762 18102 58211   100     e
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 15353 15568 15810  15660 15799 21369   100
##                                                               sd3(x, wts = w)   981   987  1022   1005  1052  1183   100
##                                                                 slow_sd(x, w)  1385  1419  3682   2023  2865 80206   100
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
##  as.centsums(x1, 1)  183  186  194    187  200  308   100  b  
##  as.centsums(x1, 2)  115  117  124    119  127  282   100 a   
##  as.centsums(x1, 3)  808  811  847    823  878  988   100   c 
##  as.centsums(x1, 4) 1400 1405 1456   1431 1491 1742   100    d
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
##  c(obj1, obj2)  15 16   20     16 18 178   100   b
##  obj3 %-% obj1  11 12   15     13 14  94   100  a
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
##    join_cent_sums(rs1, rs2) 2.0 2.1  2.6    2.2 2.4  23   100   a
##  unjoin_cent_sums(rs3, rs2) 1.8 2.0  2.1    2.1 2.2   4   100   a
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
##  as.centcosums(x1, max_ord)  55 57   64     58 70 155   100   b
##             mobj3 %-% mobj1  17 19   22     20 21  43   100  a
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
        wins, robust = TRUE), running_sum(x, wins, 
        robust = FALSE), running_mean(x, wins), roll::roll_sum(xm, 
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
##                                                           expr   min    lq  mean median    uq    max neval     cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 33376 34171 37792  35910 36477 113281   100      f 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 66933 68070 70615  69716 70904 148442   100       g
##                            running_sum(x, wins, robust = TRUE)    39    45    75     48    59   2334   100 a      
##                           running_sum(x, wins, robust = FALSE)    26    29    35     32    38     56   100 a      
##                                          running_mean(x, wins)    47    52   116     54    62   2294   100 a      
##                                       roll::roll_sum(xm, wins)  1895  1943  2012   1984  2042   2708   100  bc    
##                                      roll::roll_mean(xm, wins)  2068  2090  2169   2150  2199   2657   100  bc    
##                                        roll::roll_sd(xm, wins)  5735  5828  5936   5881  5963   7834   100    d   
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   350   371   540    387   440   2805   100 ab     
##                             RollingWindow::RollingSum(x, wins)   106   129   229    149   195   2528   100 a      
##                            RollingWindow::RollingMean(x, wins)   142   156   231    173   218   2531   100 a      
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1720  1929  2009   1965  2032   4501   100  bc    
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1730  1862  2049   2008  2061   4291   100  bc    
##                                         running_skew4(x, wins)  3436  3476  3587   3521  3581   6101   100   c    
##                                         running_kurt5(x, wins)  5535  5624  5820   5681  5821   8389   100    d   
##                                         running_tstat(x, wins)   669   678   734    691   727   3000   100 ab     
##                                       running_zscored(x, wins)   663   670   700    686   714    849   100 ab     
##                                        running_sharpe(x, wins)   663   671   691    680   702    808   100 ab     
##                                    running_apx_median(x, wins) 13001 13174 13523  13291 13541  16242   100     e  
##                                      running_centered(x, wins)   549   558   584    566   586    808   100 ab     
##                                        running_scaled(x, wins)   653   663   689    674   697    995   100 ab
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# print(resdf)
library(readr)
readr::write_csv(resdf, "timings.csv")
```

