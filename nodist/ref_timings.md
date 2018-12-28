

# fromo timings against reference implementations

Let us check timings of simple sums:


```r
library(fromo)
library(microbenchmark)

print(Sys.info())
```

```
##                                               sysname                                               release 
##                                               "Linux"                                   "4.15.0-42-generic" 
##                                               version                                              nodename 
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "67eb51117904" 
##                                               machine                                                 login 
##                                              "x86_64"                                             "unknown" 
##                                                  user                                        effective_user 
##                                              "docker"                                              "docker"
```

```r
print(sessionInfo())
```

```
## R version 3.5.0 (2018-04-23)
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6662    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       bindr_0.1.1      magrittr_1.5     grDevices_3.5.0  tidyselect_0.2.5 munsell_0.5.0    colorspace_1.3-2 R6_2.2.2         rlang_0.3.0.1   
## [10] stringr_1.3.1    plyr_1.8.4       tools_3.5.0      grid_3.5.0       gtable_0.2.0     xfun_0.4         withr_2.1.2      stats_3.5.0      lazyeval_0.2.1  
## [19] assertthat_0.2.0 tibble_1.4.2     crayon_1.3.4     bindrcpp_0.2.2   formatR_1.5      purrr_0.2.5      graphics_3.5.0   glue_1.3.0       evaluate_0.10.1 
## [28] stringi_1.2.3    compiler_3.5.0   pillar_1.3.1     scales_1.0.0     pkgconfig_2.0.2
```

```r
set.seed(12345)
x <- rnorm(1e+05)
wins <- 250

ref_running_sum <- function(x, wins) {
    cx <- cumsum(x)
    cx - c(rep(0, wins), cx[1:(length(cx) - wins)])
}

ref_running_sum2 <- function(x, wins) {
    cx <- cumsum(x)
    c(cx[1:wins], cx[(wins + 1):length(cx)] - cx[1:(length(cx) - 
        wins)])
}


# check first
blah <- running_sum(x, wins) - ref_running_sum(x, wins)
print(summary(blah[4:length(blah)]))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -2.84e-14 -3.60e-15  0.00e+00 -4.00e-16  1.80e-15  3.02e-14
```

```r
# timings
checkit <- microbenchmark(sum(x), mean(x), gc(), running_sum(x, 
    wins), running_sum(x, wins, na_rm = FALSE, restart_period = 50000L), 
    running_mean(x, wins), ref_running_sum(x, wins), 
    ref_running_sum2(x, wins), cumsum(x))

print(checkit)
```

```
## Unit: microseconds
##                                                          expr   min    lq  mean median    uq   max neval cld
##                                                        sum(x)    80    82    97     87    95   254   100 a  
##                                                       mean(x)   163   177   192    184   192   356   100 a  
##                                                          gc() 42472 42938 45354  43404 43991 95977   100   c
##                                          running_sum(x, wins)  1043  1225  1300   1280  1370  1642   100 ab 
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)  1017  1207  1311   1271  1369  1922   100 ab 
##                                         running_mean(x, wins)  1029  1181  1239   1214  1266  1851   100 ab 
##                                      ref_running_sum(x, wins)   559   773  1284   1190  1244 23254   100 ab 
##                                     ref_running_sum2(x, wins)   836  1426  1765   1723  1807 10903   100  b 
##                                                     cumsum(x)    86   111   232    266   292   548   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1299877|      1.0|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1310577|      1.0|
|ref_running_sum(x, wins)                                     | 1284122|      1.0|
|ref_running_sum2(x, wins)                                    | 1765276|      1.4|

Welford standard deviation is easy to compute quickly:




```r
library(fromo)
library(microbenchmark)

print(Sys.info())
```

```
##                                               sysname                                               release 
##                                               "Linux"                                   "4.15.0-42-generic" 
##                                               version                                              nodename 
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "67eb51117904" 
##                                               machine                                                 login 
##                                              "x86_64"                                             "unknown" 
##                                                  user                                        effective_user 
##                                              "docker"                                              "docker"
```

```r
print(sessionInfo())
```

```
## R version 3.5.0 (2018-04-23)
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
## [1] stats   utils   methods base   
## 
## other attached packages:
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6662    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       highr_0.7        pillar_1.3.1     compiler_3.5.0   formatR_1.5      plyr_1.8.4       bindr_0.1.1      tools_3.5.0      grDevices_3.5.0 
## [10] evaluate_0.10.1  tibble_1.4.2     gtable_0.2.0     lattice_0.20-35  pkgconfig_2.0.2  rlang_0.3.0.1    Matrix_1.2-14    mvtnorm_1.0-8    xfun_0.4        
## [19] withr_2.1.2      stringr_1.3.1    graphics_3.5.0   grid_3.5.0       tidyselect_0.2.5 glue_1.3.0       R6_2.2.2         survival_2.42-3  multcomp_1.4-8  
## [28] TH.data_1.0-9    purrr_0.2.5      magrittr_1.5     codetools_0.2-16 MASS_7.3-50      splines_3.5.0    scales_1.0.0     assertthat_0.2.0 colorspace_1.3-2
## [37] sandwich_2.5-0   stringi_1.2.3    lazyeval_0.2.1   munsell_0.5.0    crayon_1.3.4     zoo_1.8-4
```

```r
set.seed(12345)
x <- rnorm(1e+05)
wins <- 250

# check first
blah <- running_sd(x, wins) - ref_running_sd(x, wins)
print(summary(blah[4:length(blah)]))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -5.1e-15 -7.0e-16  1.6e-15  1.8e-15  4.7e-15  8.5e-15
```

```r
# timings
checkit <- microbenchmark(sd(x), sd3(x), ref_sd(x), 
    running_sd(x, wins), running_sd(x, wins, na_rm = FALSE, 
        restart_period = 50000L), gc(), ref_running_sd(x, 
        wins), ref_running_sd_narm(x, wins), ref_running_sd_intnel(x, 
        wins), ref_running_sd_objecty(x, wins), ref_running_sd_onecheck(x, 
        wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr    min     lq   mean median     uq    max neval   cld
##                                                        sd(x)    332    389    423    419    443    718   100 a    
##                                                       sd3(x)    849    870    909    889    930   1160   100 ab   
##                                                    ref_sd(x)    715    729    764    745    774   1024   100 ab   
##                                          running_sd(x, wins)   7981   8346   9658   8659  10655  17868   100    d 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   5777   6106   7362   6424   8205  18104   100   c  
##                                                         gc() 100783 102304 105685 103949 108026 124679   100     e
##                                      ref_running_sd(x, wins)   1182   1226   1286   1252   1328   1573   100  b   
##                                 ref_running_sd_narm(x, wins)   1175   1206   1284   1233   1303   2192   100  b   
##                               ref_running_sd_intnel(x, wins)   1182   1222   1295   1250   1363   1655   100  b   
##                              ref_running_sd_objecty(x, wins)   1186   1218   1286   1252   1352   1668   100  b   
##                             ref_running_sd_onecheck(x, wins)   1180   1218   1257   1242   1270   1488   100  b
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 9658116|      7.7|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 7362001|      5.9|
|ref_running_sd(x, wins)                                     | 1285580|      1.0|
|ref_running_sd_narm(x, wins)                                | 1284038|      1.0|
|ref_running_sd_intnel(x, wins)                              | 1295279|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1286473|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1257207|      1.0|


