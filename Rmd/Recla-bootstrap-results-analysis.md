Recla analysis - determining the LRT statistics & p-values
================
Frederick Boehm
3/7/2018

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
as_tibble(read.table("../results/2018-02-27_pvl-LD_distance_light-HP_latency-chr8-sex.txt")) -> out_0722
as_tibble(read.table("../results/2018-02-27_pvl-LD_distance_light-LD_light_pct-chr8-sex.txt")) -> out_0710
as_tibble(read.table("../results/2018-02-27_pvl-LD_light_pct-HP_latency-chr8-sex.txt")) -> out_1022
```

``` r
library(qtl2pleio)
calc_lrt_tib(out_0722) # run 552
```

    ## [1] 2.765636

``` r
calc_lrt_tib(out_0710) # run 553
```

    ## [1] 0.08584859

``` r
calc_lrt_tib(out_1022) # run 551
```

    ## [1] 2.759804

Working with the CHTC results files
-----------------------------------

Our bootstrap samples were each analyzed with computing resources from the Center for High-Throughput Computing (CHTC) at UW-Madison.

``` bash
ls ../chtc/Recla-bootstrap/results-chtc/boot-run551/*.txt | wc
ls ../chtc/Recla-bootstrap/results-chtc/boot-run552/*.txt | wc
ls ../chtc/Recla-bootstrap/results-chtc/boot-run553/*.txt | wc
```

    ##     1000    1000   74890
    ##     1000    1000   74890
    ##     1000    1000   74890

We see that each directory has the required 1000 files.

### Run 551

``` r
directory <- "../chtc/Recla-bootstrap/results-chtc/boot-run551"
fns <- dir(directory)
i <- 1
out <- list()
for (fn in fns){
  full_fn <- file.path(directory, fn)
  read.table(full_fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_1022))
```

    ## [1] 121

Run 552
-------

``` r
directory <- "../chtc/Recla-bootstrap/results-chtc/boot-run552"
fns <- dir(directory)
i <- 1
out <- list()
for (fn in fns){
  full_fn <- file.path(directory, fn)
  read.table(full_fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_0722))
```

    ## [1] 106

Run 553
-------

``` r
directory <- "../chtc/Recla-bootstrap/results-chtc/boot-run553"
fns <- dir(directory)
i <- 1
out <- list()
for (fn in fns){
  full_fn <- file.path(directory, fn)
  read.table(full_fn) -> out[[i]]
  i <- i + 1
}
sum(do.call("rbind", out) > calc_lrt_tib(out_0710))
```

    ## [1] 876
