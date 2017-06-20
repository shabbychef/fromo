

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
##                                "7b725d3a7a86"                                      "x86_64"                                     "unknown" 
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
##  [7] dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.6530      
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
##       sum(x)    80    84   102     94   108   283   100 a     
##      mean(x)   163   179   218    201   250   490   100 a     
##        sd(x)   460   513   602    568   641  1429   100 a     
##  skewness(x)  8283  9278 10515  10508 11290 14006   100   cd  
##  kurtosis(x)  8139  9005 10067  10131 10912 13744   100   c   
##       sd3(x)  1115  1225  1355   1334  1439  1936   100  b    
##     skew4(x)  9110  9666 10905  10984 11810 15077   100    d  
##     kurt5(x) 15937 17003 18774  19127 20282 22624   100     e 
##     dumbk(x) 17201 19147 21356  21645 23050 26645   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 19.6 20.7 21.2   21.1 21.6 24.6   100   c
##                                                               sd3(x, wts = w)  1.2  1.3  1.4    1.4  1.6  1.8   100 a  
##                                                                 slow_sd(x, w)  1.6  1.7  2.8    2.5  3.4  6.0   100  b
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
##  as.centsums(x1, 1)   83   87  116     96  126  408   100 a   
##  as.centsums(x1, 2)  150  163  197    176  220  394   100  b  
##  as.centsums(x1, 3)  975 1093 1213   1145 1274 2270   100   c 
##  as.centsums(x1, 4) 1664 1787 2019   1920 2228 3026   100    d
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
##  c(obj1, obj2)  17 19   25     20 21 254   100   a
##  obj3 %-% obj1  13 14   20     16 17 189   100   a
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
##    join_cent_sums(rs1, rs2) 2.5 2.8  3.4    2.9 3.1  35   100   a
##  unjoin_cent_sums(rs3, rs2) 2.1 2.4  2.7    2.6 2.7  15   100   a
##  unjoin_cent_sums(rs3, rs1) 2.1 2.4  3.2    2.5 2.7  64   100   a
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
##  as.centcosums(x1, max_ord)  51 61   80     71 91 179   100   b
##             mobj3 %-% mobj1  17 21   28     23 31  76   100  a
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
##                                                           expr   min    lq  mean median    uq    max neval      cld
##                         silly_fun(x, wins, sum, na.rm = FALSE) 31611 34487 40648  37582 45379 106662   100       g 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 64195 67572 78341  72404 85113 187635   100        h
##                                           running_sum(x, wins)   106   113   156    123   140   2530   100 a       
##                                          running_mean(x, wins)   106   116   160    127   152   2348   100 a       
##                                       roll::roll_sum(xm, wins)  1903  1993  2302   2145  2442   5987   100 ab      
##                                      roll::roll_mean(xm, wins)  2079  2169  2521   2364  2718   5231   100  b      
##                                        roll::roll_sd(xm, wins)  5749  5864  6878   6293  7194  13742   100   cd    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   360   440   663    478   598   4570   100 ab      
##                             RollingWindow::RollingSum(x, wins)   129   160   315    200   258   3259   100 ab      
##                            RollingWindow::RollingMean(x, wins)   152   201   318    244   322   2914   100 ab      
##                             RollingWindow::RollingStd(x, wins)   240   286   375    325   395   2893   100 ab      
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1763  1954  2219   2096  2423   5874   100 ab      
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1736  1996  2258   2126  2455   5577   100 ab      
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 10558 11497 14499  13547 14916 116304   100     e   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   743   832  1002    911  1043   2743   100 ab      
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   797   878  1048    960  1137   2362   100 ab      
##                                           running_sd3(x, wins)  1207  1343  1697   1488  1889   4022   100 ab      
##                                          running_skew(x, wins)  4243  4425  5149   4967  5475  13869   100   c     
##                                         running_skew4(x, wins)  4498  4725  5444   5311  5995   7524   100   c     
##                                          running_kurt(x, wins)  6056  6223  7117   7031  7881   9056   100   cd    
##                                         running_kurt5(x, wins)  6828  7130  8211   8277  8976  13748   100    d    
##                                         running_tstat(x, wins)  1029  1116  1317   1198  1389   3553   100 ab      
##                                       running_zscored(x, wins)  1018  1081  1231   1147  1279   2179   100 ab      
##                                        running_sharpe(x, wins)  1011  1065  1257   1147  1297   5643   100 ab      
##                                    running_apx_median(x, wins) 14789 15125 17119  16401 18847  22184   100      f  
##                                      running_centered(x, wins)  1099  1176  1402   1309  1557   3461   100 ab      
##                                        running_scaled(x, wins)  1024  1088  1320   1188  1460   2921   100 ab
```

```r
resdf <- rbind(resdf, checkit)
```


```r
# fake it: resdf <- readr::read_csv('timings.csv')
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
    mutate(relchange = normalized/first_norm, last_status = last_norm/first_norm) %>% 
    mutate(gtype = ifelse(grepl("^(roll|RollingWindow|RcppRoll)::", 
        expr), "brand_x", ifelse(grepl("^running_", 
        expr), "running", "summarizing")))

library(ggplot2)
ph <- allt %>% filter(!grepl("brand_x", gtype)) %>% 
    ggplot(aes(sernum, normalized, group = expr, color = expr)) + 
    geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + facet_grid(gtype ~ ., scales = "free") + 
    labs(x = "release", y = "mean time taken, relative to sum(x)", 
        title = "fromo microbenchmark timings, lasagna")
```

```
## Error in grepl("brand_x", gtype): object 'gtype' not found
```

```r
print(ph)
```

```
## Error in print(ph): object 'ph' not found
```

```r
ph <- allt %>% filter(!grepl("brand_x", gtype)) %>% 
    ggplot(aes(sernum, relchange, group = expr, color = expr)) + 
    geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + facet_grid(gtype ~ ., scales = "free") + 
    labs(x = "release", y = "normalized time taken, relative to first iteration", 
        title = "fromo microbenchmark timings, spaghetti")
```

```
## Error in grepl("brand_x", gtype): object 'gtype' not found
```

```r
print(ph)
```

```
## Error in print(ph): object 'ph' not found
```

```r
ph <- allt %>% filter(!grepl("brand_x", gtype)) %>% 
    ggplot(aes(sernum, relchange)) + geom_boxplot(aes(group = sernum), 
    alpha = 0.7) + stat_summary(aes(group = "1", color = "mean"), 
    fun.y = mean, geom = "line") + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    facet_grid(gtype ~ ., scales = "free") + labs(x = "release", 
    y = "normalized time taken, relative to first iteration", 
    color = "stat", title = "fromo microbenchmark timings, boxplots")
```

```
## Error in grepl("brand_x", gtype): object 'gtype' not found
```

```r
print(ph)
```

```
## Error in print(ph): object 'ph' not found
```

```r
allt %>% filter(!grepl("brand_x", gtype)) %>% select(expr, 
    sernum, relchange, last_status) %>% tidyr::spread(key = "sernum", 
    value = "relchange") %>% arrange(desc(last_status)) %>% 
    select(-last_status) %>% head(n = 50) %>% kable()
```

```
## Error in grepl("brand_x", gtype): object 'gtype' not found
```

```r
allt %>% filter(!grepl("brand_x", gtype)) %>% select(-sernum, 
    -relchange, -meantime, -sumx_time, -normalized) %>% 
    distinct(expr, .keep_all = TRUE) %>% arrange(desc(last_status)) %>% 
    head(n = 20) %>% kable()
```

```
## Error in grepl("brand_x", gtype): object 'gtype' not found
```


