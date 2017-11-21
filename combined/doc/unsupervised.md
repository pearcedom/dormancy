Unsupervised Analysis
================
Dominic Pearce

``` r
library(tidyverse)
library(Biobase)
library(ggthemes)
library(cowplot); theme_set(theme_grey())
source("/Volumes/igmm/sims-lab/Dominic/functions/idReplace.R")
source("/Volumes/igmm/sims-lab/Dominic/functions/mostVar.R")
source("/Volumes/igmm/sims-lab/Dominic/functions/library/mdsArrange.R")
```

<!--
//TODO try using only A, B, C samples to see if they sit close to one another
-->
``` r
dormset <- read_rds("../output/dormset.rds")
pData(dormset)[dormset$patient == 343,]$is_dormant <- TRUE
pData(dormset)[dormset$patient == 374,]$is_dormant <- TRUE
```

#### as well as how well the are able to cluster patients based on dormancy status

``` r
plotMDS <- function(eset, time_point = "all", title = ""){
    if(time_point[[1]] != "all"){
               eset_in = eset[, which(eset$timepoint %in% time_point)]
           } else {
               eset_in = eset 
           }
           arng_dfr <- mdsArrange(exprs(eset_in))
           mrg_dfr <- merge(arng_dfr, pData(eset_in), by.x = 'ids', by.y = 'xpr_id')

           p_mds <- ggplot(mrg_dfr, aes(x = x, y = y, colour = is_dormant, size = 2.1)) + 
               geom_point() + 
               labs(title = title,
                    subtitle = paste0("x separation p = ", 
                                      wilcox.test(x~is_dormant, data = mrg_dfr)$p.value)) +
               theme_pander() +
               theme(legend.position = 'top')

           p_box <- ggplot(mrg_dfr, aes(x = is_dormant, y = x, fill = is_dormant)) + 
               geom_boxplot(notch = TRUE) + 
               coord_flip() + 
               theme_pander() + 
               theme(legend.position = 'none')
           plot_grid(p_mds, p_box, rel_heights = c(3, 1), ncol = 1)
}
```

#### We can calculate the 500 most-variable genes independently for long-term treated and untreated samples

``` r
getTimepoint500 <- function(timepoint){
        dormset_sub <- dormset[, which(dormset$timepoint == timepoint)]
        mv500 <- mostVar(exprs(dormset_sub), 500) %>% row.names()
        dormset[mv500,]
}

lt_500 <- getTimepoint500("long-term") 
pre_500 <- getTimepoint500("diagnosis")

common <- intersect(row.names(lt_500), row.names(pre_500))
common_327 <- dormset[common, ]

lt_common <- which(row.names(lt_500) %in% common)
lt_unique <- row.names(lt_500)[-lt_common]
ltu_173 <- dormset[lt_unique, ]

pre_common <- which(row.names(pre_500) %in% common)
pre_unique <- row.names(pre_500)[-pre_common]
preu_173 <- dormset[pre_unique, ]
```

### Long-term 500

``` r
plotMDS(lt_500, 'all', 'All Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(lt_500, c('on-treatment', 'long-term'), 'Treated Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(lt_500, 'long-term', 'Long-term Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-3.png" style="display: block; margin: auto;" />

``` r
plotMDS(lt_500, 'on-treatment', 'On-treatment Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-4.png" style="display: block; margin: auto;" />

``` r
plotMDS(lt_500, 'diagnosis', 'Diagnostic Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-5.png" style="display: block; margin: auto;" />

### Pre-treatment 500

``` r
plotMDS(pre_500, 'all', 'All Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(pre_500, c('on-treatment', 'long-term'), 'Treated Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(pre_500, 'long-term', 'Long-term Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-3.png" style="display: block; margin: auto;" />

``` r
plotMDS(pre_500, 'on-treatment', 'On-treatment Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-4.png" style="display: block; margin: auto;" />

``` r
plotMDS(pre_500, 'diagnosis', 'Diagnostic Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-5.png" style="display: block; margin: auto;" />

### Intersect 327

``` r
plotMDS(common_327, 'all', 'All Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(common_327, c('on-treatment', 'long-term'), 'Treated Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(common_327, 'long-term', 'Long-term Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-3.png" style="display: block; margin: auto;" />

``` r
plotMDS(common_327, 'on-treatment', 'On-treatment Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-4.png" style="display: block; margin: auto;" />

``` r
plotMDS(common_327, 'diagnosis', 'Diagnostic Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-5.png" style="display: block; margin: auto;" />

### Long-term unique 173

``` r
plotMDS(ltu_173, 'all', 'All Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(ltu_173, c('on-treatment', 'long-term'), 'Treated Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(ltu_173, 'long-term', 'Long-term Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-3.png" style="display: block; margin: auto;" />

``` r
plotMDS(ltu_173, 'on-treatment', 'On-treatment Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-4.png" style="display: block; margin: auto;" />

``` r
plotMDS(ltu_173, 'diagnosis', 'Diagnostic Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-5.png" style="display: block; margin: auto;" />

### Pre-treatment unique 173

``` r
plotMDS(preu_173, 'all', 'All Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(preu_173, c('on-treatment', 'long-term'), 'Treated Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(preu_173, 'long-term', 'Long-term Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-3.png" style="display: block; margin: auto;" />

``` r
plotMDS(preu_173, 'on-treatment', 'On-treatment Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-4.png" style="display: block; margin: auto;" />

``` r
plotMDS(preu_173, 'diagnosis', 'Diagnostic Samples')
```

<img src="unsupervised_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-5.png" style="display: block; margin: auto;" />
