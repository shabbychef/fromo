

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                                "zinc" 
##                                               machine                                                 login 
##                                              "x86_64"                                                "spav" 
##                                                  user                                        effective_user 
##                                                "spav"                                                "spav"
```

```r
print(sessionInfo())
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
##  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] microbenchmark_1.4-2.1 moments_0.14           dplyr_0.7.8            ggplot2_3.1.0          fromo_0.1.3.6663       knitr_1.15.1          
##  [7] colorout_1.1-2         fortunes_1.5-4         devtools_1.13.4        drat_0.1.2            
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       bindr_0.1.1      magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0    colorspace_1.3-2 R6_2.3.0         rlang_0.3.0.1    stringr_1.3.1   
## [10] plyr_1.8.4       tools_3.3.2      grid_3.3.2       gtable_0.2.0     withr_2.1.2      lazyeval_0.2.1   digest_0.6.18    assertthat_0.2.0 tibble_1.4.2    
## [19] bindrcpp_0.2.2   formatR_1.5      purrr_0.2.5      glue_1.3.0       memoise_1.1.0    evaluate_0.10    stringi_1.2.3    pillar_1.2.1     scales_1.0.0    
## [28] pkgconfig_2.0.2
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
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -2.8e-14 -3.6e-15  0.0e+00 -4.0e-16  1.8e-15  3.0e-14
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
##                                                        sum(x)    80    81    99     83    91   373   100 a  
##                                                       mean(x)   163   178   201    183   188   665   100 a  
##                                                          gc() 35404 36166 41435  36808 37733 95089   100   c
##                                          running_sum(x, wins)  1013  1062  1207   1156  1190  2472   100 ab 
##  running_sum(x, wins, na_rm = FALSE, restart_period = 50000L)  1010  1044  1544   1115  1176 36713   100 ab 
##                                         running_mean(x, wins)  1027  1064  1255   1163  1215  2346   100 ab 
##                                      ref_running_sum(x, wins)  1086  1155  1729   1579  2123  6762   100 ab 
##                                     ref_running_sum2(x, wins)  1477  2030  3025   2460  3056 39845   100  b 
##                                                     cumsum(x)   446   461   662    569   603  1872   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sum", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                         |   meant| timeover|
|:------------------------------------------------------------|-------:|--------:|
|running_sum(x, wins)                                         | 1206545|      1.0|
|running_sum(x, wins, na_rm = FALSE, restart_period = 50000L) | 1543812|      1.3|
|ref_running_sum(x, wins)                                     | 1729013|      1.4|
|ref_running_sum2(x, wins)                                    | 3025208|      2.5|

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
## "#45~16.04.1-Ubuntu SMP Mon Nov 19 13:02:27 UTC 2018"                                                "zinc" 
##                                               machine                                                 login 
##                                              "x86_64"                                                "spav" 
##                                                  user                                        effective_user 
##                                                "spav"                                                "spav"
```

```r
print(sessionInfo())
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
##  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2         microbenchmark_1.4-2.1 moments_0.14           dplyr_0.7.8            ggplot2_3.1.0          fromo_0.1.3.6663      
##  [7] knitr_1.15.1           colorout_1.1-2         fortunes_1.5-4         devtools_1.13.4        drat_0.1.2            
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0       highr_0.6        pillar_1.2.1     formatR_1.5      plyr_1.8.4       bindr_0.1.1      tools_3.3.2      digest_0.6.18    evaluate_0.10   
## [10] memoise_1.1.0    tibble_1.4.2     gtable_0.2.0     lattice_0.20-33  pkgconfig_2.0.2  rlang_0.3.0.1    Matrix_1.2-10    mvtnorm_1.0-6    withr_2.1.2     
## [19] stringr_1.3.1    grid_3.3.2       tidyselect_0.2.5 glue_1.3.0       R6_2.3.0         survival_2.41-3  multcomp_1.4-6   TH.data_1.0-8    purrr_0.2.5     
## [28] magrittr_1.5     codetools_0.2-15 MASS_7.3-47      splines_3.3.2    scales_1.0.0     assertthat_0.2.0 colorspace_1.3-2 sandwich_2.3-4   stringi_1.2.3   
## [37] lazyeval_0.2.1   munsell_0.5.0    zoo_1.8-4
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
## -5.1e-15 -7.0e-16  1.6e-15  1.8e-15  4.7e-15  8.6e-15
```

```r
# timings
checkit <- microbenchmark(sd(x), sd3(x), ref_sd(x), 
    ref_sd_objecty(x), running_sd(x, wins), running_sd(x, 
        wins, na_rm = FALSE, restart_period = 50000L), 
    gc(), ref_running_sd(x, wins), ref_running_sd_narm(x, 
        wins), ref_running_sd_intnel(x, wins), ref_running_sd_objecty(x, 
        wins), ref_running_sd_onecheck(x, wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr   min    lq  mean median    uq   max neval  cld
##                                                        sd(x)   459   493   514    504   542   628   100 a   
##                                                       sd3(x)   848   856   893    882   920  1061   100 a   
##                                                    ref_sd(x)   714   719   749    731   778   855   100 a   
##                                            ref_sd_objecty(x)   718   723   760    752   787   859   100 a   
##                                          running_sd(x, wins) 10146 10504 11178  10674 11198 15031   100   c 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  8031  8350  9680   8484  8846 92187   100  b  
##                                                         gc() 76553 78271 79570  78969 81091 83692   100    d
##                                      ref_running_sd(x, wins)  1178  1212  1304   1246  1291  3062   100 a   
##                                 ref_running_sd_narm(x, wins)  1417  1446  1534   1496  1549  3065   100 a   
##                               ref_running_sd_intnel(x, wins)  1180  1207  1366   1249  1313  2909   100 a   
##                              ref_running_sd_objecty(x, wins)  1181  1212  1356   1259  1325  2860   100 a   
##                             ref_running_sd_onecheck(x, wins)  1181  1203  1329   1235  1313  2780   100 a
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 1.1e+07|      8.6|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 9.7e+06|      7.4|
|ref_running_sd(x, wins)                                     | 1.3e+06|      1.0|
|ref_running_sd_narm(x, wins)                                | 1.5e+06|      1.2|
|ref_running_sd_intnel(x, wins)                              | 1.4e+06|      1.1|
|ref_running_sd_objecty(x, wins)                             | 1.4e+06|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1.3e+06|      1.0|


