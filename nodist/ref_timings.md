

# fromo timings against reference implementations

Let us check timings of simple sums:


```r
library(fromo)
library(microbenchmark)

print(Sys.info())
```

```
##                                        sysname                                        release                                        version 
##                                        "Linux"                            "4.4.0-137-generic" "#163-Ubuntu SMP Mon Sep 24 13:14:43 UTC 2018" 
##                                       nodename                                        machine                                          login 
##                                 "fddb3f2044c0"                                       "x86_64"                                      "unknown" 
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6668    
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
##                                                          expr   min    lq  mean median    uq   max neval  cld
##                                                        sum(x)    76    77    83     79    86   122   100 a   
##                                                       mean(x)   155   169   175    174   184   207   100 a   
##                                                          gc() 37353 37825 38614  38294 38714 45065   100    d
##                                          running_sum(x, wins)   976  1164  1221   1230  1279  1515   100  b  
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)   976  1145  1215   1209  1296  1470   100  b  
##                                         running_mean(x, wins)   993  1132  1165   1156  1217  1449   100  b  
##                                      ref_running_sum(x, wins)   520   661  1158   1136  1210 16090   100  b  
##                                     ref_running_sum2(x, wins)   758  1277  1560   1613  1676  6163   100   c 
##                                                     cumsum(x)    79    99   212    255   275   377   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1220880|      1.1|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1214884|      1.1|
|ref_running_sum(x, wins)                                     | 1158426|      1.0|
|ref_running_sum2(x, wins)                                    | 1559738|      1.4|

Welford standard deviation is easy to compute quickly:




```r
library(fromo)
library(microbenchmark)

print(Sys.info())
```

```
##                                        sysname                                        release                                        version 
##                                        "Linux"                            "4.4.0-137-generic" "#163-Ubuntu SMP Mon Sep 24 13:14:43 UTC 2018" 
##                                       nodename                                        machine                                          login 
##                                 "fddb3f2044c0"                                       "x86_64"                                      "unknown" 
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
## [1] stats   utils   methods base   
## 
## other attached packages:
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6668    
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
##                                                         expr   min    lq  mean median     uq    max neval   cld
##                                                        sd(x)   322   378   407    395    426    584   100 a    
##                                                       sd3(x)   660   673   708    703    740    838   100 a    
##                                                    ref_sd(x)   679   695   729    719    761    838   100 a    
##                                            ref_sd_objecty(x)   658   669   708    704    733    860   100 a    
##                                          running_sd(x, wins)  7183  7610  9390   7908  10941  26784   100    d 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  5548  5754  6882   5985   7609  13128   100   c  
##                                                         gc() 89728 94347 99917  99836 104283 116612   100     e
##                                      ref_running_sd(x, wins)  1119  1175  1238   1240   1276   1402   100 ab   
##                                 ref_running_sd_narm(x, wins)  1118  1164  1249   1231   1289   1857   100 ab   
##                               ref_running_sd_intnel(x, wins)  1119  1170  1234   1232   1277   1454   100 ab   
##                              ref_running_sd_objecty(x, wins)  1125  1164  1229   1216   1279   1468   100 ab   
##                             ref_running_sd_onecheck(x, wins)  1124  1170  1231   1209   1272   1497   100 ab   
##                                 ref_running_sd_fooz(x, wins)  5312  5764  7326   6208   8430  19310   100   c  
##                                 ref_running_sd_barz(x, wins)  1737  1798  1880   1869   1929   2290   100  b   
##                                 ref_running_sd_batz(x, wins)  4679  5062  6433   5383   7312  11035   100   c
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 9389861|      7.6|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 6882105|      5.6|
|ref_running_sd(x, wins)                                     | 1238433|      1.0|
|ref_running_sd_narm(x, wins)                                | 1248792|      1.0|
|ref_running_sd_intnel(x, wins)                              | 1234121|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1229361|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1230931|      1.0|
|ref_running_sd_fooz(x, wins)                                | 7325869|      6.0|
|ref_running_sd_barz(x, wins)                                | 1879961|      1.5|
|ref_running_sd_batz(x, wins)                                | 6432595|      5.2|


