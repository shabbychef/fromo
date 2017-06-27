

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
##                                "a3888b10684b"                                      "x86_64"                                     "unknown" 
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
## -8.8e-15 -3.1e-15 -4.0e-16 -6.0e-16  2.0e-15  6.7e-15
```

```r
# timings
checkit <- microbenchmark(sd(x), sd3(x), ref_sd(x), 
    running_sd(x, wins), running_sd(x, wins, na_rm = FALSE, 
        restart_period = 50000L), ref_running_sd(x, 
        wins), ref_running_sd_narm(x, wins), ref_running_sd_intnel(x, 
        wins), ref_running_sd_objecty(x, wins), ref_running_sd_onecheck(x, 
        wins))


print(checkit)
```

```
## Unit: microseconds
##                                                         expr  min   lq  mean median    uq   max neval    cld
##                                                        sd(x)  459  486   506    495   515   664   100 a     
##                                                       sd3(x)  852  868   883    871   886  1011   100 ab    
##                                                    ref_sd(x)  715  717   733    719   732   926   100 a     
##                                          running_sd(x, wins) 9289 9990 10515  10521 10738 13468   100      f
##  running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) 7217 7600  8272   8333  8505 11849   100     e 
##                                      ref_running_sd(x, wins) 1553 1567  1994   1590  1656 37625   100   cd  
##                                 ref_running_sd_narm(x, wins) 1997 2012  2434   2032  2086 38583   100    d  
##                               ref_running_sd_intnel(x, wins) 1545 1559  1593   1578  1610  1755   100  bc   
##                              ref_running_sd_objecty(x, wins) 1669 1678  1738   1709  1750  2567   100   cd  
##                             ref_running_sd_onecheck(x, wins) 1553 1568  1643   1589  1666  4705   100   c
```

```r
checkit %>% group_by(expr) %>% dplyr::summarize(meant = mean(time, 
    na.rm = TRUE)) %>% ungroup() %>% dplyr::filter(grepl("running_sd", 
    expr)) %>% mutate(timeover = meant/min(meant, na.rm = TRUE)) %>% 
    kable()
```



|expr                                                        |   meant| timeover|
|:-----------------------------------------------------------|-------:|--------:|
|running_sd(x, wins)                                         | 1.1e+07|      6.6|
|running_sd(x, wins, na_rm = FALSE, restart_period = 50000L) | 8.3e+06|      5.2|
|ref_running_sd(x, wins)                                     | 2.0e+06|      1.2|
|ref_running_sd_narm(x, wins)                                | 2.4e+06|      1.5|
|ref_running_sd_intnel(x, wins)                              | 1.6e+06|      1.0|
|ref_running_sd_objecty(x, wins)                             | 1.7e+06|      1.1|
|ref_running_sd_onecheck(x, wins)                            | 1.6e+06|      1.0|


