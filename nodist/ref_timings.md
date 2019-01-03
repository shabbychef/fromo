

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "6aa2d6858e73" 
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
## [1] microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6710    
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
##                                                        sum(x)    80    82    87     84    91   131   100 a   
##                                                       mean(x)   162   177   188    183   193   264   100 a   
##                                                          gc() 43059 43554 44073  43860 44306 48135   100    d
##                                          running_sum(x, wins)  1030  1182  1244   1239  1311  1541   100  b  
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)  1021  1172  1224   1207  1267  1562   100  b  
##                                         running_mean(x, wins)  1043  1194  1201   1207  1232  1528   100  b  
##                                      ref_running_sum(x, wins)   572   671  1163   1178  1211 15007   100  b  
##                                     ref_running_sum2(x, wins)   861  1316  1617   1690  1728  5876   100   c 
##                                                     cumsum(x)    85   109   217    264   276   365   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1243662|      1.1|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1224404|      1.1|
|ref_running_sum(x, wins)                                     | 1162545|      1.0|
|ref_running_sum2(x, wins)                                    | 1616819|      1.4|

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                        "6aa2d6858e73" 
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
## [1] bindrcpp_0.2.2       microbenchmark_1.4-6 moments_0.14         dplyr_0.7.8          ggplot2_3.1.0        knitr_1.21           fromo_0.1.3.6710    
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
## -4.4e-15 -8.0e-16  1.6e-15  1.8e-15  4.8e-15  8.5e-15
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
##                                                         expr    min     lq   mean median     uq    max neval    cld
##                                                        sd(x)    332    366    396    399    409    531   100 a     
##                                                       sd3(x)    695    699    728    703    729    919   100  b    
##                                                    ref_sd(x)    715    719    742    724    744    899   100  b    
##                                            ref_sd_objecty(x)    693    696    725    711    735    902   100  b    
##                                          running_sd(x, wins)   2951   3009   3123   3112   3185   3601   100     e 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)   1187   1254   1375   1344   1376   5441   100   c   
##                                                         gc() 102845 103633 104328 104049 104494 108825   100      f
##                                      ref_running_sd(x, wins)   1185   1231   1310   1331   1359   1628   100   c   
##                                 ref_running_sd_narm(x, wins)   1184   1222   1304   1330   1354   1545   100   c   
##                               ref_running_sd_intnel(x, wins)   1184   1223   1304   1329   1351   1680   100   c   
##                              ref_running_sd_objecty(x, wins)   1181   1214   1310   1332   1358   1836   100   c   
##                             ref_running_sd_onecheck(x, wins)   1185   1230   1313   1332   1368   1550   100   c   
##                                 ref_running_sd_fooz(x, wins)   1201   1223   1315   1342   1365   1563   100   c   
##                                 ref_running_sd_barz(x, wins)   1814   1914   1982   1986   2037   2311   100    d
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 3122525|      2.4|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 1374994|      1.1|
|ref_running_sd(x, wins)                                     | 1310291|      1.0|
|ref_running_sd_narm(x, wins)                                | 1304341|      1.0|
|ref_running_sd_intnel(x, wins)                              | 1304179|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1309712|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1313480|      1.0|
|ref_running_sd_fooz(x, wins)                                | 1314768|      1.0|
|ref_running_sd_barz(x, wins)                                | 1982323|      1.5|


