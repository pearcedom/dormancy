Timepoint correlations
================
Dominic Pearce

``` r
library(tidyverse)
library(Biobase)
library(ggforce)
```

``` r
dormset <- read_rds("../output/dormset.rds")
```

#### Here we're going to characterise our patients based on the sampling timepoints - *diagnostic*, *on-treatment* and *long-term* - where we would anticipate desensitised patients to exhibit decreasing correlation at later timepoints, and vice versa for dormant patients.

#### We can also extend this categorical comparison to simply compare correlation vs. time on treatment.

``` r
corArrange <- function(eset){
    cor_mtx <- cor(exprs(eset))
    diag(cor_mtx) <- NA
    cor_vec <- data.frame(cor = colMeans(cor_mtx, na.rm = TRUE))
    cor_mrg <- merge(cor_vec, pData(eset), by = 0)
    cor_mrg$timepoint <- factor(cor_mrg$timepoint, levels = c("diagnosis", "on-treatment", "long-term"))
    cor_mrg
}

corByStatus <- function(eset){
    lapply(c(TRUE, FALSE), function(logical){
               corArrange(dormset[, which(dormset$is_dormant == logical)])
}) %>% do.call(rbind, .)
}

all_cor <- corByStatus(dormset)
```

#### Intra-status correlations

``` r
ggplot(all_cor, aes(x = is_dormant, y = cor)) + 
    geom_boxplot(outlier.size = 0, notch = TRUE) + 
    geom_sina()
```

<img src="timepoint-correlations_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#### Intra-status correlations over time

``` r
ggplot(all_cor, aes(x = days_treated, y = cor, colour = is_dormant)) + 
    geom_point() + 
    geom_smooth()
```

<img src="timepoint-correlations_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### Intra-status & intra-timepoint correlations

``` r
timepoint_cor <- lapply(c("diagnosis", "on-treatment", "long-term"), function(timepoint){
                            corByStatus(dormset[, which(dormset$timepoint == timepoint)])
}) %>% do.call(rbind, .)

ggplot(timepoint_cor, aes(x = timepoint, y = cor)) + 
    geom_boxplot(outlier.size = 0, notch = TRUE) +
#    geom_jitter(width = 0.2) +
    geom_sina() +
    facet_wrap(~is_dormant, nrow = 1)
```

<img src="timepoint-correlations_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

#### Intra-status & intra-timepoint correlations

``` r
ggplot(timepoint_cor, aes(x = days_treated, y = cor, colour = is_dormant)) +
    geom_point() +
    geom_smooth()
```

<img src="timepoint-correlations_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

#### Check that dormancy status is correct...

``` r
pheno <- pData(dormset)
library(reshape2)
xpr <- melt(exprs(dormset))

mrg <- merge(xpr, pheno, by.x = 'Var2', by.y = 'sample_id')

mrg$timepoint <- factor(mrg$timepoint, levels = c("diagnosis", "on-treatment", "long-term"))

ggplot(mrg[which(mrg$Var1 %in% c("MKI67", "MCM2", "PCNA")),], aes(x = timepoint, y = value)) +
    geom_boxplot() +
    facet_grid(Var1~is_dormant, scales = 'free')
```

<img src="timepoint-correlations_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
