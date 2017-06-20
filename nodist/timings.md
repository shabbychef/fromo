

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
##                                "e6471c5e0c0c"                                      "x86_64"                                     "unknown" 
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
##       sum(x)    80    80    94     85    95   244   100 a     
##      mean(x)   162   165   198    181   213   500   100 a     
##        sd(x)   458   471   544    514   586   839   100 a     
##  skewness(x)  8272  8634  9816   9436 10226 15445   100   cd  
##  kurtosis(x)  8118  8223  9320   8955 10015 14288   100   c   
##       sd3(x)  1115  1129  1246   1203  1346  1676   100  b    
##     skew4(x)  9043  9238 10322   9939 10851 16662   100    d  
##     kurt5(x) 15756 15975 17834  16793 19138 26183   100     e 
##     dumbk(x) 17033 17500 19555  18912 20639 28067   100      f
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
##  cent_moments(x, max_order = 4, wts = w, na_rm = TRUE, normalize_wts = FALSE) 20.8 21.5 21.9   21.9 22.2 23.3   100   c
##                                                               sd3(x, wts = w)  1.2  1.4  1.5    1.5  1.6  1.8   100 a  
##                                                                 slow_sd(x, w)  1.6  2.1  2.7    2.5  3.3  5.3   100  b
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
##  as.centsums(x1, 1)   85   93  119    105  138  272   100 a   
##  as.centsums(x1, 2)  153  174  204    192  220  361   100  b  
##  as.centsums(x1, 3) 1025 1157 1247   1224 1304 1614   100   c 
##  as.centsums(x1, 4) 1821 1964 2101   2103 2195 2474   100    d
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
##  c(obj1, obj2)  16 18   26     19 23 302   100   a
##  obj3 %-% obj1  12 14   20     15 17 229   100   a
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
##    join_cent_sums(rs1, rs2) 2.7 3.0  4.2    3.1 3.5  68   100   a
##  unjoin_cent_sums(rs3, rs2) 2.5 2.7  3.6    2.9 3.2  39   100   a
##  unjoin_cent_sums(rs3, rs1) 2.5 2.6  3.8    2.9 3.2  79   100   a
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
##                        expr min lq mean median  uq max neval cld
##  as.centcosums(x1, max_ord)  62 70   92     80 105 302   100   b
##             mobj3 %-% mobj1  20 23   32     25  36 101   100  a
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
##                         silly_fun(x, wins, sum, na.rm = FALSE) 34317 37469 42103  41897 44157 136457   100          j 
##                        silly_fun(x, wins, mean, na.rm = FALSE) 68046 79487 83537  83065 84862 178732   100           k
##                                           running_sum(x, wins)   111   122   149    135   163    386   100 a          
##                                          running_mean(x, wins)   110   121   142    133   154    251   100 a          
##                                       roll::roll_sum(xm, wins)  1928  2057  2257   2163  2377   3200   100 a  d       
##                                      roll::roll_mean(xm, wins)  2069  2304  2571   2491  2716   5404   100    d       
##                                        roll::roll_sd(xm, wins)  5773  6161  6833   6620  7367   9119   100     e g    
##       RollingWindow::RollingSum(x, wins, na_method = "ignore")   384   446   679    532   635   3535   100 a  d       
##                             RollingWindow::RollingSum(x, wins)   109   161   256    202   241   3080   100 ab         
##                            RollingWindow::RollingMean(x, wins)   147   193   348    247   304   2731   100 a c        
##                             RollingWindow::RollingStd(x, wins)   248   300   481    357   425   2882   100 a  d       
##   RcppRoll::roll_sum(xm, n = wins, align = "right", fill = NA)  1855  2059  2342   2200  2410   6299   100  bcd       
##  RcppRoll::roll_mean(xm, n = wins, align = "right", fill = NA)  1845  2121  2375   2304  2508   4949   100   cd       
##    RcppRoll::roll_sd(xm, n = wins, align = "right", fill = NA) 11180 12653 14766  13922 15483  92290   100        h   
##    running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   780   923  1041    991  1172   1511   100 a  d       
##      running_sd(x, wins, na_rm = TRUE, restart_period = 1000L)   834   954  1093   1062  1176   1677   100 a  d       
##                                           running_sd3(x, wins)  1257  1530  1822   1654  1949   8468   100 a  d       
##                                          running_skew(x, wins)  4302  4816  5299   5031  5659  10633   100     e      
##                                         running_skew4(x, wins)  4526  5213  5746   5617  6160   8855   100     ef     
##                                          running_kurt(x, wins)  6421  6821  7486   7376  8018   9701   100      fg    
##                                         running_kurt5(x, wins)  7112  7820  8518   8517  9144  10071   100       g    
##                                         running_tstat(x, wins)  1037  1220  1388   1358  1491   2396   100 a  d       
##                                       running_zscored(x, wins)  1044  1183  1323   1256  1459   1875   100 a  d       
##                                        running_sharpe(x, wins)  1053  1184  1359   1276  1466   3188   100 a  d       
##                                    running_apx_median(x, wins) 14887 16785 18284  18726 19222  21555   100         i  
##                                      running_centered(x, wins)  1113  1273  1436   1346  1603   2415   100 a  d       
##                                        running_scaled(x, wins)  1027  1192  1368   1342  1518   1996   100 a  d
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
    mutate(gtype = case_when(grepl("^(roll|RollingWindow|RcppRoll)::", 
        .$expr) ~ "brand_x", grepl("^running_", .$expr) ~ 
        "running", TRUE ~ "summarizing"))

library(ggplot2)
ph <- allt %>% filter(!grepl("brand_x", gtype)) %>% 
    ggplot(aes(sernum, normalized, group = expr, color = expr)) + 
    geom_line() + geom_point() + scale_y_log10() + 
    guides(colour = FALSE) + theme(axis.text.x = element_text(angle = -90, 
    hjust = 0)) + facet_grid(gtype ~ ., scales = "free") + 
    labs(x = "release", y = "mean time taken, relative to sum(x)", 
        title = "fromo microbenchmark timings")
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
        title = "fromo microbenchmark timings")
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
    color = "stat", title = "fromo microbenchmark timings")
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


