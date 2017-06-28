

# fromo timings against reference implementations

Welford standard deviation is easy to compute quickly:




```r
library(fromo)
library(microbenchmark)

print(Sys.info())
```

```
##                                       sysname                                       release                                       version 
##                                       "Linux"                            "4.4.0-77-generic" "#98-Ubuntu SMP Wed Apr 26 08:34:02 UTC 2017" 
##                                      nodename                                       machine                                         login 
##                                "aa8d2266bc72"                                      "x86_64"                                     "unknown" 
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
## [1] microbenchmark_1.4-2.1 moments_0.14           dplyr_0.5.0            ggplot2_2.2.1          knitr_1.15.1           fromo_0.1.3.6560      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11     grDevices_3.3.0  assertthat_0.1   R6_2.1.2         grid_3.3.0       plyr_1.8.4       DBI_0.6-1        gtable_0.2.0     formatR_1.4     
## [10] magrittr_1.5     evaluate_0.10    scales_0.4.1     stringi_1.0-1    lazyeval_0.2.0   graphics_3.3.0   tools_3.3.0      stringr_1.0.0    munsell_0.4.3   
## [19] stats_3.3.0      colorspace_1.3-2 tibble_1.2
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
    running_sd(x, wins), running_sd(x, wins, na_rm = FALSE, 
        restart_period = 50000L), gc(), ref_running_sd(x, 
        wins), ref_running_sd_narm(x, wins), ref_running_sd_intnel(x, 
        wins), ref_running_sd_objecty(x, wins), ref_running_sd_onecheck(x, 
        wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr   min    lq  mean median    uq   max neval    cld
##                                                        sd(x)   465   492   519    504   535   888   100 a     
##                                                       sd3(x)   850   856   899    874   910  1428   100 ab    
##                                                    ref_sd(x)   715   720   760    739   789  1048   100 ab    
##                                          running_sd(x, wins)  9248 10004 10931  10567 11051 48333   100     e 
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L)  7134  7907  8304   8365  8724 10514   100    d  
##                                                         gc() 36000 36983 38062  38343 38855 42031   100      f
##                                      ref_running_sd(x, wins)  1186  1206  1259   1240  1297  1458   100  bc   
##                                 ref_running_sd_narm(x, wins)  1421  1444  1505   1475  1548  1870   100   c   
##                               ref_running_sd_intnel(x, wins)  1184  1196  1258   1218  1298  1765   100  bc   
##                              ref_running_sd_objecty(x, wins)  1182  1202  1254   1222  1277  1732   100  bc   
##                             ref_running_sd_onecheck(x, wins)  1186  1203  1257   1223  1294  1536   100  bc
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 1.1e+07|      8.7|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 8.3e+06|      6.6|
|ref_running_sd(x, wins)                                     | 1.3e+06|      1.0|
|ref_running_sd_narm(x, wins)                                | 1.5e+06|      1.2|
|ref_running_sd_intnel(x, wins)                              | 1.3e+06|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1.3e+06|      1.0|
|ref_running_sd_onecheck(x, wins)                            | 1.3e+06|      1.0|


