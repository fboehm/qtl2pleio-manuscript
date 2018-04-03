Recla analysis - determining the LRT statistics & p-values
================
Frederick Boehm
2018-04-03 18:00:12

In this document, we use the 1000 files (for each of the bivariate traits) that each contain a likelihood ratio test statistic for a bootstrap sample.

Load the scan\_pvl results for each of the 3 bivariate traits
-------------------------------------------------------------

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
as_tibble(read.table("recla-07-22.txt")) -> out_0722
as_tibble(read.table("recla-07-10.txt")) -> out_0710
as_tibble(read.table("recla-10-22.txt")) -> out_1022
```

``` r
library(qtl2pleio)
calc_lrt_tib(out_0722) # run 562
```

    ## [1] 2.759916

``` r
calc_lrt_tib(out_0710) # run 563
```

    ## [1] 0.09387654

``` r
calc_lrt_tib(out_1022) # run 561
```

    ## [1] 2.771408

Working with the CHTC results files
-----------------------------------

Our bootstrap samples were each analyzed with computing resources from the Center for High-Throughput Computing (CHTC) at UW-Madison.

``` bash
ls ../chtc/Recla-bootstrap/submit_files/*run561*.txt | wc
ls ../chtc/Recla-bootstrap/submit_files/*run562*.txt | wc
ls ../chtc/Recla-bootstrap/submit_files/*run563*.txt | wc
```

    ##     1000    1000   62890
    ##     1000    1000   62890
    ##     1000    1000   62890

We see that each analysis has the required 1000 files.

### Run 561

``` r
directory <- "../chtc/Recla-bootstrap/submit_files"
fns <- dir(directory, pattern = "recla-boot-run561", full.names = TRUE)
i <- 1
out <- list()
for (fn in fns){
  read.table(fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_1022))
```

    ## [1] 109

Run 562
-------

``` r
fns <- dir(directory, pattern = "recla-boot-run562", full.names = TRUE)
i <- 1
out <- list()
for (fn in fns){
  read.table(fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_0722))
```

    ## [1] 108

Run 563
-------

``` r
fns <- dir(directory, pattern = "recla-boot-run563", full.names = TRUE)
i <- 1
out <- list()
for (fn in fns){
  read.table(fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_0710))
```

    ## [1] 871

``` r
devtools::session_info()
```

    ## Session info -------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.3 (2017-11-30)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2018-04-03

    ## Packages -----------------------------------------------------------------

    ##  package    * version    date       source                             
    ##  assertthat   0.2.0      2017-04-11 CRAN (R 3.4.0)                     
    ##  backports    1.1.2      2017-12-13 CRAN (R 3.4.3)                     
    ##  base       * 3.4.3      2017-12-07 local                              
    ##  bindr        0.1.1      2018-03-13 CRAN (R 3.4.4)                     
    ##  bindrcpp     0.2        2017-06-17 CRAN (R 3.4.0)                     
    ##  compiler     3.4.3      2017-12-07 local                              
    ##  datasets   * 3.4.3      2017-12-07 local                              
    ##  devtools     1.13.5     2018-02-18 CRAN (R 3.4.3)                     
    ##  digest       0.6.15     2018-01-28 cran (@0.6.15)                     
    ##  dplyr      * 0.7.4      2017-09-28 cran (@0.7.4)                      
    ##  evaluate     0.10.1     2017-06-24 CRAN (R 3.4.0)                     
    ##  glue         1.2.0      2017-10-29 CRAN (R 3.4.2)                     
    ##  graphics   * 3.4.3      2017-12-07 local                              
    ##  grDevices  * 3.4.3      2017-12-07 local                              
    ##  htmltools    0.3.6      2017-04-28 CRAN (R 3.4.0)                     
    ##  knitr        1.20       2018-02-20 CRAN (R 3.4.3)                     
    ##  lubridate    1.7.3      2018-02-27 CRAN (R 3.4.3)                     
    ##  magrittr     1.5.0      2017-09-23 Github (tidyverse/magrittr@0a76de2)
    ##  memoise      1.1.0      2017-04-21 CRAN (R 3.4.0)                     
    ##  methods    * 3.4.3      2017-12-07 local                              
    ##  pillar       1.2.1      2018-02-27 CRAN (R 3.4.3)                     
    ##  pkgconfig    2.0.1      2017-03-21 CRAN (R 3.4.0)                     
    ##  qtl2pleio  * 0.1.2      2018-03-14 local (fboehm/qtl2pleio@NA)        
    ##  R6           2.2.2      2017-06-17 CRAN (R 3.4.0)                     
    ##  Rcpp         0.12.16    2018-03-13 cran (@0.12.16)                    
    ##  rlang        0.2.0.9000 2018-03-13 Github (tidyverse/rlang@49e9389)   
    ##  rmarkdown    1.9        2018-03-01 CRAN (R 3.4.3)                     
    ##  rprojroot    1.3-2      2018-01-03 CRAN (R 3.4.3)                     
    ##  stats      * 3.4.3      2017-12-07 local                              
    ##  stringi      1.1.7      2018-03-12 CRAN (R 3.4.4)                     
    ##  stringr      1.3.0      2018-02-19 cran (@1.3.0)                      
    ##  tibble       1.4.2      2018-01-22 cran (@1.4.2)                      
    ##  tools        3.4.3      2017-12-07 local                              
    ##  utils      * 3.4.3      2017-12-07 local                              
    ##  withr        2.1.2      2018-03-15 CRAN (R 3.4.4)                     
    ##  yaml         2.1.18     2018-03-08 CRAN (R 3.4.4)
