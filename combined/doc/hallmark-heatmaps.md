Hallmarks heatmaps
================
Dominic Pearce

<!--
####Setup
-->
``` r
library(tidyverse)
library(Biobase)
library(GSEABase)
library(cowplot); theme_set(theme_grey())
source("../../lib/arrangeGSEAGenesets.R")
source("../../../../functions/library/heatmapArrange.R")
```

``` r
dormset <- read_rds("../output/dormset.rds")
hallmarks <- getGmt("../data/h.all.v6.1.symbols.gmt") %>% arrangeGenesets(dormset)

#Significant hallmarks for various comparisons
hm_all <- readLines("../output/hm-gsea-unpaired-all.txt")
hm_diagnosis <- readLines("../output/hm-gsea-unpaired-diagnosis.txt")
hm_on <- readLines("../output/hm-gsea-unpaired-on-treatment.txt")
hm_long <- readLines("../output/hm-gsea-unpaired-long-term.txt")
```

#### Here we're looking to take a geneset that we've determined as informative in *hallmarks-gsea.Rmd* and display the per-sample expression as a heatmap.

``` r
heatmapPlot <- function(eset, geneset, title){
        #arrange
        esub <- eset[geneset,]
        dfr <- heatmapArrange(exprs(esub), cluster_row = TRUE, cluster_column = FALSE)
        mrg_dfr <- merge(dfr, pData(esub), by.x = 'variable', by.y = 'xpr_id')
        ord_levels <- unique(mrg_dfr$variable[order(mrg_dfr$days_treated)])
        mrg_dfr$ordered_samples <- factor(mrg_dfr$variable, levels = ord_levels)
        fill_max <- range(mrg_dfr$value) %>% abs() %>% max()

        #plot
        #p_hmap <- 
                ggplot(mrg_dfr, aes(x = ordered_samples, y = row_value, fill = value)) +
                        geom_tile() + 
                        #scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850') +
                        scale_fill_gradientn(values=c(1, .7, .5, .3, 0), 
                                             colours=c("#fb6a4a", 
                                                           "#cb181d", 
                                                           "black", 
                                                           "#238b45", 
                                                           "#74c476"), 
                                             limits = c(-fill_max, fill_max),
                                             na.value = "white") +

                        facet_grid(is_dormant~timepoint, scales = 'free_x', space = 'free') + 
                        labs(title = title) +
                        theme(
                              axis.text.x = element_blank(),
                              axis.ticks.x = element_blank()
                              )

        #p_bar <- ggplot(mrg_dfr, aes(x = variable, y = "Status", fill = is_dormant)) +
        #                geom_tile() +
        #                theme(
        #                      axis.text.x = element_blank(),
        #                      axis.ticks.x = element_blank()
        #                      ) + 
        #                labs(title = title)
        #plot_grid(p_bar, p_hmap, ncol = 1, align = 'v', rel_heights = c(1, 8))
}
```

#### All sample comparisons

``` r
lapply(hm_all, function(x) heatmapPlot(dormset, hallmarks[[x]], x))
```

    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"

    ## [[1]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

    ## 
    ## [[2]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-2.png" style="display: block; margin: auto;" />

    ## 
    ## [[3]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-3.png" style="display: block; margin: auto;" />

    ## 
    ## [[4]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-4.png" style="display: block; margin: auto;" />

    ## 
    ## [[5]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-5.png" style="display: block; margin: auto;" />

    ## 
    ## [[6]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-6.png" style="display: block; margin: auto;" />

    ## 
    ## [[7]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-7.png" style="display: block; margin: auto;" />

    ## 
    ## [[8]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-8.png" style="display: block; margin: auto;" />

#### Diagnostic comparisons

``` r
lapply(hm_diagnosis, function(x) {
               eset <- dormset[, which(dormset$timepoint == "diagnosis")]
               heatmapPlot(eset, hallmarks[[x]], x)
      })
```

    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"

    ## [[1]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

    ## 
    ## [[2]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-2.png" style="display: block; margin: auto;" />

    ## 
    ## [[3]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-3.png" style="display: block; margin: auto;" />

    ## 
    ## [[4]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-4.png" style="display: block; margin: auto;" />

    ## 
    ## [[5]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-5.png" style="display: block; margin: auto;" />

    ## 
    ## [[6]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-6.png" style="display: block; margin: auto;" />

    ## 
    ## [[7]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-7.png" style="display: block; margin: auto;" />

    ## 
    ## [[8]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-8.png" style="display: block; margin: auto;" />

    ## 
    ## [[9]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-9.png" style="display: block; margin: auto;" />

    ## 
    ## [[10]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-10.png" style="display: block; margin: auto;" />

    ## 
    ## [[11]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-11.png" style="display: block; margin: auto;" />

    ## 
    ## [[12]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-12.png" style="display: block; margin: auto;" />

    ## 
    ## [[13]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-13.png" style="display: block; margin: auto;" />

    ## 
    ## [[14]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-14.png" style="display: block; margin: auto;" />

    ## 
    ## [[15]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-15.png" style="display: block; margin: auto;" />

#### On-treatment comparisons

``` r
lapply(hm_on, function(x) {
               eset <- dormset[, which(dormset$timepoint == "on-treatment")]
               heatmapPlot(eset, hallmarks[[x]], x)
      })
```

    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"

    ## [[1]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

    ## 
    ## [[2]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-2.png" style="display: block; margin: auto;" />

    ## 
    ## [[3]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-3.png" style="display: block; margin: auto;" />

#### Long-term comparisons

``` r
lapply(hm_long, function(x) {
               eset <- dormset[, which(dormset$timepoint == "long-term")]
               heatmapPlot(eset, hallmarks[[x]], x)
      })
```

    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"
    ## [1] "scale_fill_gradient2 reminder : scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')"

    ## [[1]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

    ## 
    ## [[2]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

    ## 
    ## [[3]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-3.png" style="display: block; margin: auto;" />

    ## 
    ## [[4]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-4.png" style="display: block; margin: auto;" />

    ## 
    ## [[5]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-5.png" style="display: block; margin: auto;" />

    ## 
    ## [[6]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-6.png" style="display: block; margin: auto;" />

    ## 
    ## [[7]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-7.png" style="display: block; margin: auto;" />

    ## 
    ## [[8]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-8.png" style="display: block; margin: auto;" />

    ## 
    ## [[9]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-9.png" style="display: block; margin: auto;" />

    ## 
    ## [[10]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-10.png" style="display: block; margin: auto;" />

    ## 
    ## [[11]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-11.png" style="display: block; margin: auto;" />

    ## 
    ## [[12]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-12.png" style="display: block; margin: auto;" />

    ## 
    ## [[13]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-13.png" style="display: block; margin: auto;" />

    ## 
    ## [[14]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-14.png" style="display: block; margin: auto;" />

    ## 
    ## [[15]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-15.png" style="display: block; margin: auto;" />

    ## 
    ## [[16]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-16.png" style="display: block; margin: auto;" />

    ## 
    ## [[17]]

<img src="hallmark-heatmaps_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-17.png" style="display: block; margin: auto;" />
