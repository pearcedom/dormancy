Batch effects and correction
================
Dominic Pearce

``` r
library(tidyverse)
library(Biobase)
library(sva)
library(reshape2)
library(ggthemes)
library(cowplot); theme_set(theme_gray())
source("../../../../functions/mostVar.R")
source("../../../../functions/library/mdsArrange.R")
```

#### We have our frma and loess normalised expression set, but there still clearly exists a big batch difference between FF and FFPE source material. This is evident in terms of the overall sample distributions

``` r
georgeset_no_loess <- read_rds("../output/georgeset.rds")
georgeset <- read_rds("../output/georgeset-frma.rds")
```

``` r
batchEffectDists <- function(mtx_input){
    #arrange
    george_mlt <- melt(mtx_input)
    george_mrg <- merge(george_mlt, pData(georgeset), by.x = "Var2", by.y = 'xpr_id')
    george_mrg$Var2 <- factor(george_mrg$Var2, levels = unique(george_mrg$Var2[order(george_mrg$is_freshfrozen)]))
    #boxplot
    p_box <- ggplot(george_mrg, aes(x = Var2, y = value, fill = is_freshfrozen)) + 
        geom_boxplot() +
        theme(legend.position = 'bottom',
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    #densityplot
    p_dist <- ggplot(george_mrg, aes(x = value, colour = is_freshfrozen, group = Var2)) + 
        geom_density() + 
        theme_pander() +    
        theme(legend.position = 'none',
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    #combine
    plot_grid(p_box, p_dist, ncol = 1)
}

batchEffectDists(exprs(georgeset))
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

#### ...as well as when comparing just the 500 most variable genes as an mds

``` r
batchEffectMDS <- function(mtx_input, colour_by){
    #arrange
    mv500 <- mostVar(mtx_input, 500) %>% row.names()
    arg500 <- mdsArrange(mtx_input[mv500,]) 
    mds_input <- merge(arg500, pData(georgeset), by.x = 'ids', by.y = 'xpr_id')
    #and plot
    ggplot(mds_input, aes_string("x", "y", colour = colour_by)) + 
        geom_point() + 
        theme_pander() + 
        theme(legend.position = 'bottom')
}

batchEffectMDS(exprs(georgeset), "is_freshfrozen")
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#### We can additionally check a number of other variables to check for any other effects

``` r
var_vec <- c("biopsy", "is_uss", "Age", "ER", "Diagnostic.CB.Grade", "No.of.Pos.Nodes")
plot_lst <- lapply(var_vec, function(variable) batchEffectMDS(exprs(georgeset), variable))
cowplot::plot_grid(plotlist = plot_lst, ncol = 3)
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### And as it's perhaps the most likely source of biological variation, we'll plot biopsy number a little larger

``` r
batchEffectMDS(exprs(georgeset), "timepoint")
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

#### We'll attempt to correct for this source material-specific batch effect using ComBat

``` r
g_raw <- exprs(georgeset)
g_cb <- ComBat(exprs(georgeset), batch = georgeset$is_freshfrozen)
```

    ## Standardizing Data across genes

``` r
#source("../../../../functions/xpn-conor.R")
#g_cb <- xpn(exprs(georgeset)[georgeset$is_freshfrozen], exprs(georgeset)[!georgeset$is_freshfrozen])
```

``` r
batchEffectDists(g_cb)
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

#### Better? What about by MDS?

``` r
batchEffectMDS(g_cb, "is_freshfrozen")
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

#### Still a clear divide by source material but let's check the FF-FFPE replicate samples only

``` r
material_rep_ids <- c("50-1", "50-2", "50-3", "50-4", "135-1", "135-4", "188-1", "268-1", "268-2", "298-1", "298-2", "347-3", "413-1", "413-2", "416-1", "416-2")
material_rep_vec <- rep(FALSE, ncol(georgeset))
material_rep_vec[which(georgeset$sample_id %in% material_rep_ids)] <- TRUE

batchEffectDists(g_cb[, material_rep_vec])
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(g_cb[, material_rep_vec], "is_freshfrozen")
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

#### Supplementary

For reference the non-loess normalised, rma normalised data

``` r
batchEffectDists(exprs(georgeset_no_loess))
```

<img src="/Volumes/igmm/sims-lab/Dominic/labbook/dormancy/georgetown/doc/normalisation-batch-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />
