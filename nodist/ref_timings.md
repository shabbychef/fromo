

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "cf8e3d5f6221" 
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6703    
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
##                                                          expr   min    lq  mean median    uq    max neval cld
##                                                        sum(x)    80    82    99     86    94    374   100  a 
##                                                       mean(x)   162   178   209    183   201    575   100  a 
##                                                          gc() 41598 42097 51791  42795 45829 101659   100   b
##                                          running_sum(x, wins)  1031  1235  1334   1309  1382   2301   100  a 
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)  1021  1216  1338   1328  1404   2009   100  a 
##                                         running_mean(x, wins)  1046  1193  1291   1226  1311   2023   100  a 
##                                      ref_running_sum(x, wins)   558   761  1306   1188  1285  15262   100  a 
##                                     ref_running_sum2(x, wins)   836  1472  1993   1699  1921   6269   100  a 
##                                                     cumsum(x)    85   116   249    265   285    967   100  a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1333799|      1.0|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1338434|      1.0|
|ref_running_sum(x, wins)                                     | 1305601|      1.0|
|ref_running_sum2(x, wins)                                    | 1993476|      1.5|

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "cf8e3d5f6221" 
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
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6703    
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
        wins), ref_running_sd_barz(x, wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr   min    lq   mean median     uq    max neval    cld
##                                                        sd(x)   330   361    383    388    395    481   100 a     
##                                                       sd3(x)   694   699    725    714    744    812   100 ab    
##                                                    ref_sd(x)   715   717    739    724    743    899   100  b    
##                                            ref_sd_objecty(x)   692   695    712    697    715    870   100 ab    
##                                          running_sd(x, wins)  3564  3650   3718   3711   3772   4107   100     e 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  1801  1884   1981   1949   1991   5835   100    d  
##                                                         gc() 98476 99269 100625  99697 100812 117881   100      f
##                                      ref_running_sd(x, wins)  1180  1199   1290   1330   1353   1499   100   c   
##                                 ref_running_sd_narm(x, wins)  1179  1202   1279   1270   1333   1583   100   c   
##                               ref_running_sd_intnel(x, wins)  1181  1206   1278   1259   1333   1525   100   c   
##                              ref_running_sd_objecty(x, wins)  1178  1203   1290   1326   1356   1574   100   c   
##                             ref_running_sd_onecheck(x, wins)  1180  1210   1298   1328   1351   1513   100   c   
##                                 ref_running_sd_fooz(x, wins)  1806  1837   1924   1923   1989   2451   100    d  
##                                 ref_running_sd_barz(x, wins)  1803  1900   1971   1972   2034   2248   100    d
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 3718146|      2.9|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 1980530|      1.6|
|ref_running_sd(x, wins)                                     | 1290156|      1.0|
|ref_running_sd_narm(x, wins)                                | 1279373|      1.0|
|ref_running_sd_intnel(x, wins)                              | 1278020|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1289746|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1298289|      1.0|
|ref_running_sd_fooz(x, wins)                                | 1923960|      1.5|
|ref_running_sd_barz(x, wins)                                | 1971476|      1.5|


