

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "c304cdc5a241" 
##                                               machine                                                 login 
##                                              "x86_64"                                             "unknown" 
##                                                  user                                        effective_user 
##                                              "docker"                                              "docker"
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6700    
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
##                                                        sum(x)    80    82    90     84    92   145   100 a  
##                                                       mean(x)   162   176   185    181   189   323   100 a  
##                                                          gc() 41709 42073 43990  42660 43763 59568   100   c
##                                          running_sum(x, wins)  1058  1216  1261   1262  1314  1837   100  b 
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)  1018  1204  1260   1257  1332  1625   100  b 
##                                         running_mean(x, wins)  1045  1193  1209   1205  1232  1526   100  b 
##                                      ref_running_sum(x, wins)   557   736  1208   1184  1229 16764   100  b 
##                                     ref_running_sum2(x, wins)   862  1304  1646   1691  1738  7386   100  b 
##                                                     cumsum(x)    85   110   221    263   278   419   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1260878|      1.0|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1259554|      1.0|
|ref_running_sum(x, wins)                                     | 1207673|      1.0|
|ref_running_sum2(x, wins)                                    | 1645951|      1.4|

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "c304cdc5a241" 
##                                               machine                                                 login 
##                                              "x86_64"                                             "unknown" 
##                                                  user                                        effective_user 
##                                              "docker"                                              "docker"
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
## [1] stats   utils   methods base   
## 
## other attached packages:
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6700    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       highr_0.7        pillar_1.3.1     compiler_3.5.2   formatR_1.5      plyr_1.8.4       bindr_0.1.1      tools_3.5.2      grDevices_3.5.2 
## [10] evaluate_0.12    tibble_1.4.2     gtable_0.2.0     lattice_0.20-38  pkgconfig_2.0.2  rlang_0.3.0.1    Matrix_1.2-15    mvtnorm_1.0-8    xfun_0.4        
## [19] withr_2.1.2      stringr_1.3.1    graphics_3.5.2   grid_3.5.2       tidyselect_0.2.5 glue_1.3.0       R6_2.3.0         survival_2.43-3  multcomp_1.4-8  
## [28] TH.data_1.0-9    purrr_0.2.5      magrittr_1.5     codetools_0.2-15 MASS_7.3-51.1    splines_3.5.2    scales_1.0.0     assertthat_0.2.0 colorspace_1.3-2
## [37] sandwich_2.5-0   stringi_1.2.4    lazyeval_0.2.1   munsell_0.5.0    crayon_1.3.4     zoo_1.8-4
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
## -5.1e-15 -7.0e-16  1.6e-15  1.8e-15  4.7e-15  8.8e-15
```

```r
# timings
checkit <- microbenchmark(sd(x), sd3(x), ref_sd(x), 
    ref_sd_objecty(x), running_sd(x, wins), running_sd(x, 
        wins, na_rm = FALSE, restart_period = 50000L), 
    gc(), ref_running_sd(x, wins), ref_running_sd_narm(x, 
        wins), ref_running_sd_intnel(x, wins), ref_running_sd_objecty(x, 
        wins), ref_running_sd_onecheck(x, wins), ref_running_sd_fooz(x, 
        wins), ref_running_sd_barz(x, wins), ref_running_sd_batz(x, 
        wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr   min     lq   mean median     uq    max neval cld
##                                                        sd(x)   331    383    429    408    441    837   100 a  
##                                                       sd3(x)   695    707    753    727    757   1070   100 a  
##                                                    ref_sd(x)   715    727    771    742    769   1122   100 a  
##                                            ref_sd_objecty(x)   693    696    751    717    753   1077   100 a  
##                                          running_sd(x, wins)  3442   3529   3794   3605   3804   5991   100 ab 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  1661   1700   1824   1750   1864   3155   100 a  
##                                                         gc() 99139 100409 116553 108055 114500 236523   100   c
##                                      ref_running_sd(x, wins)  1184   1214   1297   1253   1342   1803   100 a  
##                                 ref_running_sd_narm(x, wins)  1179   1207   1348   1231   1346   2573   100 a  
##                               ref_running_sd_intnel(x, wins)  1182   1219   1324   1251   1354   1972   100 a  
##                              ref_running_sd_objecty(x, wins)  1186   1213   1287   1242   1311   1878   100 a  
##                             ref_running_sd_onecheck(x, wins)  1180   1213   1312   1253   1343   1921   100 a  
##                                 ref_running_sd_fooz(x, wins)  1668   1713   1891   1756   1873   3368   100 a  
##                                 ref_running_sd_barz(x, wins)  1841   1880   2007   1936   2021   3001   100 a  
##                                 ref_running_sd_batz(x, wins)  5437   5613   7198   5940   7479  19362   100  b
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 3793716|      3.0|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 1823973|      1.4|
|ref_running_sd(x, wins)                                     | 1297370|      1.0|
|ref_running_sd_narm(x, wins)                                | 1347750|      1.1|
|ref_running_sd_intnel(x, wins)                              | 1324354|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1286822|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1311987|      1.0|
|ref_running_sd_fooz(x, wins)                                | 1890950|      1.5|
|ref_running_sd_barz(x, wins)                                | 2007416|      1.6|
|ref_running_sd_batz(x, wins)                                | 7198232|      5.6|


