Recla analysis: Updated figures
================
Frederick Boehm
4/16/2018

Read pvl scan results from files
--------------------------------

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
as_tibble(read.table("recla-07-10.txt")) -> pvl0710
as_tibble(read.table("recla-07-22.txt")) -> pvl0722
as_tibble(read.table("recla-10-22.txt")) -> pvl1022
```

Load Recla from qtl2data
------------------------

``` r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
recla <- read_cross2(file)
# make sex a covariate for use in pvl_scan
recla[[6]][ , 1, drop = FALSE] -> sex
# insert pseudomarkers
insert_pseudomarkers(recla, step = 0.10) -> pseudomap
gm <- pseudomap$`8`
```

``` r
probs <- calc_genoprob(recla, map = pseudomap)
```

We now convert the genotype probabilities to haplotype dosages.

``` r
aprobs <- genoprob_to_alleleprob(probs)
```

We now calculate kinship matrices, by the "leave one chromosome out (loco)" method.

``` r
kinship <- calc_kinship(aprobs, "loco")
```

``` r
recla$pheno -> ph
log(ph) -> lph
apply(FUN = broman::winsorize, X = lph, MARGIN = 2) -> wlph
```

``` r
sex2 <- sex == "female"
out <- scan1(genoprobs = aprobs, pheno = wlph, kinship = kinship, addcovar = sex2, reml = TRUE)
```

``` r
(peaks <- find_peaks(out, pseudomap, threshold = 5) %>%
  arrange(chr, pos) %>%
   select(- lodindex))
```

    ##                    lodcolumn chr     pos      lod
    ## 1                         bw   1 23.7675 5.246220
    ## 2         OF_distance_first4   1 43.2385 5.795778
    ## 3                OF_distance   2 49.9770 5.396383
    ## 4            OF_immobile_pct   2 51.6240 9.542721
    ## 5                         bw   2 52.3932 7.371267
    ## 6  VC_bottom_distance_first4   2 71.0160 6.536997
    ## 7         OF_distance_first4   3 10.7360 5.587546
    ## 8  VC_bottom_distance_first4   3 16.3700 5.523711
    ## 9            VC_top_time_pct   3 17.9390 5.955033
    ## 10         LD_distance_light   3 23.4390 5.251984
    ## 11                        bw   3 24.8390 5.698274
    ## 12        VC_top_time_first4   3 48.1280 6.053878
    ## 13           VC_top_velocity   3 48.5630 6.257610
    ## 14             OF_corner_pct   4  8.9111 6.436684
    ## 15           OF_immobile_pct   4 37.5206 5.507398
    ## 16         LD_distance_light   4 71.2992 5.025103
    ## 17                        bw   5 10.0740 5.785036
    ## 18 VC_bottom_distance_first4   5 19.9741 5.405890
    ## 19        VC_top_time_first4   5 20.0741 6.032911
    ## 20        VC_bottom_distance   5 20.5930 5.649634
    ## 21        VC_bottom_time_pct   5 20.5930 6.337760
    ## 22       TS_latency_immobile   5 43.3504 6.028778
    ## 23             OF_corner_pct   5 64.2551 5.665595
    ## 24           OF_immobile_pct   6 53.4292 6.967008
    ## 25     TS_frequency_climbing   6 57.0362 5.267422
    ## 26                        bw   7  9.1778 6.036994
    ## 27          OF_periphery_pct   7 41.7912 5.706906
    ## 28          TS_time_immobile   7 49.6778 8.038799
    ## 29        OF_distance_first4   7 57.9454 5.411283
    ## 30           VC_top_distance   7 83.8778 5.713105
    ## 31           OF_immobile_pct   7 83.9778 5.808357
    ## 32     TS_frequency_climbing   8 48.1732 5.487263
    ## 33         LD_distance_light   8 55.2762 5.322643
    ## 34              LD_light_pct   8 55.2762 5.268480
    ## 35                HP_latency   8 57.7732 6.205672
    ## 36         LD_distance_light   9 36.6965 5.143105
    ## 37              LD_light_pct   9 36.6965 5.392595
    ## 38        VC_top_time_first4   9 38.4834 5.098868
    ## 39           VC_top_time_pct   9 39.2680 6.397931
    ## 40                HP_latency   9 46.8502 5.228065
    ## 41                        bw  10  3.7781 6.500065
    ## 42        OF_distance_first4  10 29.6698 5.457855
    ## 43        VC_bottom_time_pct  10 32.5438 5.458701
    ## 44                HP_latency  10 60.6000 5.110472
    ## 45          OF_periphery_pct  10 74.8530 5.205981
    ## 46           VC_top_distance  11  7.8200 6.252911
    ## 47    VC_top_distance_first4  11 11.6236 5.181317
    ## 48 VC_bottom_distance_first4  11 54.3420 5.393215
    ## 49            LD_transitions  11 60.1040 5.624933
    ## 50     VC_bottom_transitions  11 60.5984 5.319996
    ## 51              LD_light_pct  11 63.3943 6.465399
    ## 52         LD_distance_light  11 63.4514 6.394151
    ## 53           VC_top_time_pct  12 20.6776 6.961159
    ## 54        VC_bottom_velocity  12 21.7760 5.655678
    ## 55             OF_center_pct  12 35.5140 6.399420
    ## 56          OF_periphery_pct  12 53.5776 7.228409
    ## 57             OF_corner_pct  13 59.7966 6.594594
    ## 58     VC_bottom_time_first4  14 11.9183 5.090795
    ## 59 VC_bottom_distance_first4  14 12.5316 5.592832
    ## 60     VC_bottom_transitions  14 12.5316 6.325557
    ## 61           VC_top_velocity  14 12.7819 6.822648
    ## 62        VC_bottom_distance  14 14.5316 5.377262
    ## 63     TS_frequency_climbing  14 21.1141 5.424521
    ## 64             OF_center_pct  14 53.7316 5.377182
    ## 65     TS_frequency_climbing  15 12.6680 6.043675
    ## 66              LD_light_pct  15 15.2374 5.643566
    ## 67        OF_distance_first4  16 23.2656 5.245592
    ## 68           VC_top_distance  17 15.6390 6.669809
    ## 69           VC_top_velocity  18  8.3750 5.608030
    ## 70        VC_top_time_first4  18 18.1850 6.477911
    ## 71            LD_transitions  18 37.4182 5.069268
    ## 72 VC_bottom_distance_first4  19 24.9615 7.369096
    ## 73     VC_bottom_time_first4  19 24.9615 7.479033
    ## 74           OF_immobile_pct  19 31.9505 5.595291
    ## 75                HP_latency  19 47.7977 5.402501

``` r
peaks8 <- peaks %>%
  filter(chr == 8, pos > 50, pos < 60
         )
pos_LD_distance_light <- peaks8 %>%
  filter(lodcolumn == "LD_distance_light") %>%
  select(pos)
pos_LD_light_pct <- peaks8 %>%
  filter(lodcolumn == "LD_light_pct") %>%
  select(pos)
pos_HP_latency <- peaks8 %>%
  filter(lodcolumn == "HP_latency") %>%
  select(pos)
```

``` r
library(xtable)
print(xtable(peaks), include.rownames = FALSE, tabular.environment = "longtable")
```

    ## Warning in print.xtable(xtable(peaks), include.rownames = FALSE,
    ## tabular.environment = "longtable"): Attempt to use "longtable" with
    ## floating = TRUE. Changing to FALSE.

    ## % latex table generated in R 3.4.3 by xtable 1.8-2 package
    ## % Tue Apr 17 13:25:24 2018
    ## \begin{longtable}{llrr}
    ##   \hline
    ## lodcolumn & chr & pos & lod \\ 
    ##   \hline
    ## bw & 1 & 23.77 & 5.25 \\ 
    ##   OF\_distance\_first4 & 1 & 43.24 & 5.80 \\ 
    ##   OF\_distance & 2 & 49.98 & 5.40 \\ 
    ##   OF\_immobile\_pct & 2 & 51.62 & 9.54 \\ 
    ##   bw & 2 & 52.39 & 7.37 \\ 
    ##   VC\_bottom\_distance\_first4 & 2 & 71.02 & 6.54 \\ 
    ##   OF\_distance\_first4 & 3 & 10.74 & 5.59 \\ 
    ##   VC\_bottom\_distance\_first4 & 3 & 16.37 & 5.52 \\ 
    ##   VC\_top\_time\_pct & 3 & 17.94 & 5.96 \\ 
    ##   LD\_distance\_light & 3 & 23.44 & 5.25 \\ 
    ##   bw & 3 & 24.84 & 5.70 \\ 
    ##   VC\_top\_time\_first4 & 3 & 48.13 & 6.05 \\ 
    ##   VC\_top\_velocity & 3 & 48.56 & 6.26 \\ 
    ##   OF\_corner\_pct & 4 & 8.91 & 6.44 \\ 
    ##   OF\_immobile\_pct & 4 & 37.52 & 5.51 \\ 
    ##   LD\_distance\_light & 4 & 71.30 & 5.03 \\ 
    ##   bw & 5 & 10.07 & 5.79 \\ 
    ##   VC\_bottom\_distance\_first4 & 5 & 19.97 & 5.41 \\ 
    ##   VC\_top\_time\_first4 & 5 & 20.07 & 6.03 \\ 
    ##   VC\_bottom\_distance & 5 & 20.59 & 5.65 \\ 
    ##   VC\_bottom\_time\_pct & 5 & 20.59 & 6.34 \\ 
    ##   TS\_latency\_immobile & 5 & 43.35 & 6.03 \\ 
    ##   OF\_corner\_pct & 5 & 64.26 & 5.67 \\ 
    ##   OF\_immobile\_pct & 6 & 53.43 & 6.97 \\ 
    ##   TS\_frequency\_climbing & 6 & 57.04 & 5.27 \\ 
    ##   bw & 7 & 9.18 & 6.04 \\ 
    ##   OF\_periphery\_pct & 7 & 41.79 & 5.71 \\ 
    ##   TS\_time\_immobile & 7 & 49.68 & 8.04 \\ 
    ##   OF\_distance\_first4 & 7 & 57.95 & 5.41 \\ 
    ##   VC\_top\_distance & 7 & 83.88 & 5.71 \\ 
    ##   OF\_immobile\_pct & 7 & 83.98 & 5.81 \\ 
    ##   TS\_frequency\_climbing & 8 & 48.17 & 5.49 \\ 
    ##   LD\_distance\_light & 8 & 55.28 & 5.32 \\ 
    ##   LD\_light\_pct & 8 & 55.28 & 5.27 \\ 
    ##   HP\_latency & 8 & 57.77 & 6.21 \\ 
    ##   LD\_distance\_light & 9 & 36.70 & 5.14 \\ 
    ##   LD\_light\_pct & 9 & 36.70 & 5.39 \\ 
    ##   VC\_top\_time\_first4 & 9 & 38.48 & 5.10 \\ 
    ##   VC\_top\_time\_pct & 9 & 39.27 & 6.40 \\ 
    ##   HP\_latency & 9 & 46.85 & 5.23 \\ 
    ##   bw & 10 & 3.78 & 6.50 \\ 
    ##   OF\_distance\_first4 & 10 & 29.67 & 5.46 \\ 
    ##   VC\_bottom\_time\_pct & 10 & 32.54 & 5.46 \\ 
    ##   HP\_latency & 10 & 60.60 & 5.11 \\ 
    ##   OF\_periphery\_pct & 10 & 74.85 & 5.21 \\ 
    ##   VC\_top\_distance & 11 & 7.82 & 6.25 \\ 
    ##   VC\_top\_distance\_first4 & 11 & 11.62 & 5.18 \\ 
    ##   VC\_bottom\_distance\_first4 & 11 & 54.34 & 5.39 \\ 
    ##   LD\_transitions & 11 & 60.10 & 5.62 \\ 
    ##   VC\_bottom\_transitions & 11 & 60.60 & 5.32 \\ 
    ##   LD\_light\_pct & 11 & 63.39 & 6.47 \\ 
    ##   LD\_distance\_light & 11 & 63.45 & 6.39 \\ 
    ##   VC\_top\_time\_pct & 12 & 20.68 & 6.96 \\ 
    ##   VC\_bottom\_velocity & 12 & 21.78 & 5.66 \\ 
    ##   OF\_center\_pct & 12 & 35.51 & 6.40 \\ 
    ##   OF\_periphery\_pct & 12 & 53.58 & 7.23 \\ 
    ##   OF\_corner\_pct & 13 & 59.80 & 6.59 \\ 
    ##   VC\_bottom\_time\_first4 & 14 & 11.92 & 5.09 \\ 
    ##   VC\_bottom\_distance\_first4 & 14 & 12.53 & 5.59 \\ 
    ##   VC\_bottom\_transitions & 14 & 12.53 & 6.33 \\ 
    ##   VC\_top\_velocity & 14 & 12.78 & 6.82 \\ 
    ##   VC\_bottom\_distance & 14 & 14.53 & 5.38 \\ 
    ##   TS\_frequency\_climbing & 14 & 21.11 & 5.42 \\ 
    ##   OF\_center\_pct & 14 & 53.73 & 5.38 \\ 
    ##   TS\_frequency\_climbing & 15 & 12.67 & 6.04 \\ 
    ##   LD\_light\_pct & 15 & 15.24 & 5.64 \\ 
    ##   OF\_distance\_first4 & 16 & 23.27 & 5.25 \\ 
    ##   VC\_top\_distance & 17 & 15.64 & 6.67 \\ 
    ##   VC\_top\_velocity & 18 & 8.38 & 5.61 \\ 
    ##   VC\_top\_time\_first4 & 18 & 18.18 & 6.48 \\ 
    ##   LD\_transitions & 18 & 37.42 & 5.07 \\ 
    ##   VC\_bottom\_distance\_first4 & 19 & 24.96 & 7.37 \\ 
    ##   VC\_bottom\_time\_first4 & 19 & 24.96 & 7.48 \\ 
    ##   OF\_immobile\_pct & 19 & 31.95 & 5.60 \\ 
    ##   HP\_latency & 19 & 47.80 & 5.40 \\ 
    ##    \hline
    ## \hline
    ## \end{longtable}

Plots
-----

``` r
library(qtl2pleio)
p1022 <- tidy_scan_pvl(pvl1022, pmap = gm) %>%
  add_intercepts(c(as.numeric(pos_LD_light_pct), as.numeric(pos_HP_latency))) %>%
  plot_pvl(phenames = colnames(recla$pheno)[c(10, 22)])
```

    ## Warning: Column `marker1`/`marker` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `marker2`/`marker` joining factor and character vector,
    ## coercing into character vector

``` r
p0722 <- tidy_scan_pvl(pvl0722, pmap = gm) %>%
  add_intercepts(c(as.numeric(pos_LD_distance_light), as.numeric(pos_HP_latency))) %>%
  plot_pvl(phenames = colnames(recla$pheno)[c(7, 22)])
```

    ## Warning: Column `marker1`/`marker` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `marker2`/`marker` joining factor and character vector,
    ## coercing into character vector

``` r
p0710 <- tidy_scan_pvl(pvl0710, pmap = gm) %>%
  add_intercepts(c(as.numeric(pos_LD_distance_light), as.numeric(pos_LD_light_pct))) %>%
  plot_pvl(phenames = colnames(recla$pheno)[c(7, 10)])
```

    ## Warning: Column `marker1`/`marker` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `marker2`/`marker` joining factor and character vector,
    ## coercing into character vector

``` r
library(cowplot)
```

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'ggplot2'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     vars

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
myplot <- plot_grid(p1022, p0722, p0710, labels=c("A", "B", "C"), ncol = 1)
```

    ## Warning: Removed 208 rows containing missing values (geom_path).

    ## Warning: Removed 209 rows containing missing values (geom_path).

``` r
save_plot(filename = "all3.png", plot = myplot, base_aspect_ratio = 2.5, ncol = 1, nrow = 3)
```

Three scatter plots of phenotypes
---------------------------------

``` r
as_tibble(wlph) -> wlph_tib
scatter0710 <- ggplot() + geom_point(data = wlph_tib, aes(y = LD_distance_light, x = LD_light_pct))
scatter1022 <- ggplot() + geom_point(data = wlph_tib, aes(y = HP_latency, x = LD_light_pct))
scatter0722 <- ggplot() + geom_point(data = wlph_tib, aes(y = HP_latency, x = LD_distance_light))
myplot_scatter <- plot_grid(scatter1022, scatter0722, scatter0710, labels=c("A", "B", "C"), ncol = 1)
```

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

``` r
save_plot(filename = "all3scatter.png", plot = myplot_scatter, base_aspect_ratio = 2.5, ncol = 1, nrow = 3)
```

Genome-wide LOD plots for the 3 traits from Recla
-------------------------------------------------

``` r
png("genomewide_lod_trait7.png")
plot(out, map = pseudomap, lodcolumn = 7)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("genomewide_lod_trait10.png")
plot(out, map = pseudomap, lodcolumn = 10)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("genomewide_lod_trait22.png")
plot(out, map = pseudomap, lodcolumn = 22)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Allele effects plots on chr 8 for each of the three Recla traits
----------------------------------------------------------------

``` r
sex <- 
png("coef7.png")
scan1coef(aprobs[ , 8], pheno = wlph[, 7], kinship = kinship$`8`, 
          reml = TRUE,
          addcovar = sex2) %>%
  plot_coefCC(map = pseudomap, scan1_output = out[, 7, drop = FALSE],
              top_panel_prop = 0.5)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("coef10.png")
scan1coef(aprobs[ , 8], pheno = wlph[, 10], kinship = kinship$`8`, 
          reml = TRUE,
          addcovar = sex2) %>%
  plot_coefCC(map = pseudomap, scan1_output = out[, 10, drop = FALSE],
              top_panel_prop = 0.5)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("coef22.png")
scan1coef(aprobs[ , 8], pheno = wlph[, 22], kinship = kinship$`8`, 
          reml = TRUE,
          addcovar = sex2) %>%
  plot_coefCC(map = pseudomap, scan1_output = out[, 22, drop = FALSE],
              top_panel_prop = 0.5)
dev.off()
```

    ## quartz_off_screen 
    ##                 2
