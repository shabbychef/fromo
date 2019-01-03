

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
##                                 "149ad6718c93"                                       "x86_64"                                      "unknown" 
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6701    
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
##                                                        sum(x)    76    78    85     80    86   160   100 a  
##                                                       mean(x)   154   168   178    175   186   271   100 a  
##                                                          gc() 37407 38014 39520  38612 40286 56851   100   c
##                                          running_sum(x, wins)   988  1165  1227   1232  1297  1466   100  b 
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)   979  1151  1225   1218  1290  1665   100  b 
##                                         running_mean(x, wins)  1000  1134  1168   1150  1224  1388   100  b 
##                                      ref_running_sum(x, wins)   530   676  1149   1131  1215 14440   100  b 
##                                     ref_running_sum2(x, wins)   840  1335  1623   1627  1733  5727   100  b 
##                                                     cumsum(x)    79   105   217    258   274   424   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1226529|      1.1|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1224693|      1.1|
|ref_running_sum(x, wins)                                     | 1149266|      1.0|
|ref_running_sum2(x, wins)                                    | 1622627|      1.4|

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
##                                 "149ad6718c93"                                       "x86_64"                                      "unknown" 
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
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6701    
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
##                                                         expr   min    lq  mean median    uq    max neval      cld
##                                                        sd(x)   314   349   372    367   394    495   100 a       
##                                                       sd3(x)   659   664   693    680   723    787   100 a       
##                                                    ref_sd(x)   679   682   713    700   733    808   100 ab      
##                                            ref_sd_objecty(x)   658   660   688    679   707    748   100 a       
##                                          running_sd(x, wins)  3247  3316  3402   3389  3468   3768   100      f  
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  1573  1601  1666   1664  1711   1884   100    de   
##                                                         gc() 89451 90044 92014  90441 92196 108353   100        h
##                                      ref_running_sd(x, wins)  1117  1134  1188   1173  1225   1361   100  bc     
##                                 ref_running_sd_narm(x, wins)  1118  1136  1182   1162  1197   1818   100  bc     
##                               ref_running_sd_intnel(x, wins)  1118  1137  1178   1174  1199   1402   100  bc     
##                              ref_running_sd_objecty(x, wins)  1116  1134  1181   1162  1217   1426   100  bc     
##                             ref_running_sd_onecheck(x, wins)  1119  1143  1192   1185  1233   1400   100   cd    
##                                 ref_running_sd_fooz(x, wins)  1586  1604  1676   1675  1723   1926   100     e   
##                                 ref_running_sd_barz(x, wins)  1717  1785  1838   1834  1874   2132   100     e   
##                                 ref_running_sd_batz(x, wins)  4739  4921  5671   5060  5347  11232   100       g
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 3402379|      2.9|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 1665753|      1.4|
|ref_running_sd(x, wins)                                     | 1187544|      1.0|
|ref_running_sd_narm(x, wins)                                | 1182196|      1.0|
|ref_running_sd_intnel(x, wins)                              | 1177541|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1181283|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1192342|      1.0|
|ref_running_sd_fooz(x, wins)                                | 1676443|      1.4|
|ref_running_sd_barz(x, wins)                                | 1837674|      1.6|
|ref_running_sd_batz(x, wins)                                | 5670503|      4.8|


