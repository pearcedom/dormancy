Hallmarks GSEA
================
Dominic Pearce

<!--
//TODO  -When having two scores with the same name its average was used. 
        -Zeros were removed. 
        -Scores without names (which can not be in any signature) removed. 
        -Non complete cases (i.e. NAs, NaNs) were removed
-->
``` r
library(tidyverse)
library(GSEABase)
library(phenoTest)
library(testthat)
library(ggthemes)
```

``` r
dormset <- read_rds("../output/dormset.rds")
```

#### Here we're going to compare enrichment of the [Hallmark genesets](https://www.ncbi.nlm.nih.gov/pubmed/26771021) between dormant and desensitised patients.

### Functions

#### First we'll write a fold-change function for expressionsets

``` r
getFC <- function(eset, id1, id2, id_col, patient_col){
    intset <- eset[, which(pData(eset)[[id_col]] == id1 | pData(eset)[[id_col]] == id2)]
    singles <- table(intset[[patient_col]]) %>% .[.==1] %>% names()
    fcset <- intset[, -which(intset[[patient_col]] %in% singles)]
    a_set <- fcset[, fcset[[id_col]] == id1]
    b_set <- fcset[, fcset[[id_col]] == id2]
    exprs(b_set) - exprs(a_set)
}
```

#### GSEA calculation functions (for unpaired and paired comparisons)

``` r
arrangeGSEA <- function(gsea_out, comparison = NA, dormancy = NA){
    dfr <- data.frame(summary(gsea_out)[, c("nes", "fdr")])
    dfr[dfr$fdr >= 0.05, "nes"] <- NA
    dfr$hallmark <- row.names(dfr)
    dfr$comparison <- comparison
    dfr$dormancy <- dormancy
    dfr$direction <- ifelse(grepl("UP$", dfr$hallmark), "UP", 
           ifelse(grepl("DN$", dfr$hallmark), "DOWN", "EITHER"))
    dfr
}

arrangePairedGSEA <- function(geneset_list){
    #for domants
    eset_base <- dormset[, which(dormset$is_dormant)]
    eset_dorm <- eset_base[,  -which(eset_base$timepoint_split != "A")]

    dorm_dlt <- getFC(eset_dorm, "diagnosis", "long-term", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "D-LT", dormancy = "Dormant")

    dorm_dot <- getFC(eset_dorm, "diagnosis", "on-treatment", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "D-OT", dormancy = "Dormant")

    dorm_otlt <- getFC(eset_dorm, "on-treatment", "long-term", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "OT-LT", dormancy = "Dormant")

    #for desensitiseds
    eset_base2 <- dormset[, which(!dormset$is_dormant)]
    eset_dssn <- eset_base2[,  -which(eset_base2$timepoint_split != "A")]

    dssn_dlt <- getFC(eset_dssn, "diagnosis", "long-term", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "D-LT", dormancy = "Desensitised")

    dssn_dot <- getFC(eset_dssn, "diagnosis", "on-treatment", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "D-OT", dormancy = "Desensitised")

    dssn_otlt <- getFC(eset_dssn, "on-treatment", "long-term", "timepoint", "patient") %>% 
                rowMeans() %>% 
                sort() %>% 
                gsea(gsets = geneset_list, logScale = FALSE, center = TRUE) %>% 
                arrangeGSEA(comparison = "OT-LT", dormancy = "Desensitised")

    #arrange for ggplot
    gsea_dfr <- do.call(rbind, list(dorm_dlt, dorm_dot, dorm_otlt, dssn_dlt, dssn_dot, dssn_otlt))
    gsea_dfr$dormancy <- factor(gsea_dfr$dormancy, levels = c("Dormant", "Desensitised")) 
    gsea_dfr$comparison <- factor(gsea_dfr$comparison, levels = c("D-OT", "D-LT", "OT-LT")) 
    gsea_dfr
}
```

#### Geneset arrangement

``` r
arrangeGenesets <- function(GMT){
    genesets <- geneIds(GMT)
    intersets <- lapply(genesets, function(set){
                            intersect(set, row.names(dormset))
    })
    intersets
}
```

### Analysis

#### First we'll assemble the Hallmark and cancer-specific genesets, as a sensible easy to use list of ids that are present in `dormset`.

``` r
hallmarks <- getGmt("../data/h.all.v6.1.symbols.gmt") %>% arrangeGenesets()
cancersets <- getGmt("../data/c6.all.v6.1.symbols.gmt") %>% arrangeGenesets()
```

#### We'll begin with an unpaired comparison, by arranging our data as two sets - dormants & desensitiseds - calculating the significant **NES** scores, and plotting.

``` r
arrangeUnpairedGSEA <- function(eset, geneset){
    #separate into dormants and desensitiseds
    dorm_xpr <- exprs(eset[, which(eset$is_dormant)]) %>% rowMeans()
    dssd_xpr <- exprs(eset[, which(!eset$is_dormant)]) %>% rowMeans()
    #calculate fold-change - this orientation is relative to desensitiseds
    fc_vec <- sort(dssd_xpr - dorm_xpr)

    test_that("test that our vector meets the phenoTest reccomended specs", {
        expect_false(unique(fc_vec %>% is.na()))
        expect_false(unique(fc_vec == 0))
        #even at 15 decimals points, for all genes there are value ties but they're such a small
        #percentage of all genes that I don't think it is critical to the results here
    })
    #perform GSEA and arrange for plotting
    geneset_out <- gsea(x = fc_vec, gsets = geneset, logScale = FALSE)
    arrangeGSEA(geneset_out, comparison = "Dormant-Desensitised")
} 


plotUnpaired <- function(x){
    ggplot(x, aes(x = comparison, y = hallmark, fill = nes)) + 
            geom_tile(colour = "white") +
            scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white") +
            labs(title = "GSEA : Dormant vs. Desensitised", subtitle = "Enrichment is relative to acquired resistant patients") 
}
```

``` r
arrangeUnpairedGSEA(dormset, hallmarks) %>%
    plotUnpaired()
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
set.seed(125)
arrangeUnpairedGSEA(dormset, cancersets) %>%
        .[!is.na(.$nes),] %>%
        plotUnpaired() +
        scale_fill_gradient2(high = "#fe9929", mid = "white", low = "#4eb3d3", na.value = "white") + 
        facet_grid(direction~., scales = 'free', space = 'free')
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

#### We can then split our dormant vs. desensitised comparison by timepoint

``` r
lapply(c("diagnosis", "on-treatment", "long-term"), function(x) {
    dfr <- arrangeUnpairedGSEA(dormset[, which(dormset$timepoint == x)], hallmarks) 
    dfr$timepoint <- factor(x, levels = c("diagnosis", "on-treatment", "long-term"))
    dfr
    }) %>% 
    do.call(rbind, .) %>%
    ggplot(aes(x = timepoint, y = hallmark, fill = nes)) +
        geom_tile(colour = "#f0f0f0") + 
        scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white") +
        #facet_wrap(~timepoint, nrow = 1) +
        labs(title = "Hallmark genesets") + 
        theme_pander()
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
lapply(c("diagnosis", "on-treatment", "long-term"), function(x) {
    dfr <- arrangeUnpairedGSEA(dormset[, which(dormset$timepoint == x)], cancersets) 
    dfr$timepoint <- factor(x, levels = c("diagnosis", "on-treatment", "long-term"))
    dfr
    }) %>% 
    do.call(rbind, .) %>%
    ggplot(aes(x = timepoint, y = hallmark, fill = nes)) +
        geom_tile(colour = "#f0f0f0") + 
        scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white") +
        #facet_wrap(~timepoint, nrow = 1) +
        labs(title = "Hallmark genesets") + 
#        theme_pander() +
        facet_grid(direction~., scales = 'free', space = 'free')
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

#### Now we'll take advantage of our paired data and look at fold change overtime, separtely for dormant and desensitised patients

#### Combine and plot

**you were finishing getting back to the same place with hm 50 (need to run the function defined above and then plot) and then going to apply the C6 gensets**

**also add clustering to the function above?**

``` r
hm_paired <- arrangePairedGSEA(hallmarks)
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

``` r
ggplot(hm_paired, aes(x = dormancy, y = hallmark, fill = nes)) +
    geom_tile(colour = "#f0f0f0") +
    facet_wrap(~comparison, nrow = 1) +
    scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white")
```

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

``` r
ggplot(hm_paired, aes(x = comparison, y = hallmark, fill = nes)) +
    geom_tile(colour = "#f0f0f0") +
    facet_wrap(~dormancy, nrow = 1) +
    scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white")
```

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-2.png" style="display: block; margin: auto;" />

#### Enrichment differences between Dormant and Desensitised in &gt;= 2 comparisons

-   xenobiotic metabolism.
-   uv response up
-   wnt signalling
-   ros
-   interferon alpha

``` r
cs_paired <- arrangePairedGSEA(cancersets)
```

    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.
    ## The provided gene sets have more than 10 distinct lengths, therefore gam approximation will be used.

``` r
ggplot(cs_paired, aes(x = comparison, y = hallmark, fill = nes)) +
    geom_tile(colour = "#f0f0f0") +
    facet_wrap(~dormancy, nrow = 1) +
    scale_fill_gradientn(values=c(1, .7, .5, .3, 0), colours=rev(c("#4eb3d3", "#a8ddb5", "white", "#fee391", "#fe9929")), na.value = "white") +
    facet_grid(direction~., scales = 'free', space = 'free')
```

<img src="hallmark-gsea_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />
