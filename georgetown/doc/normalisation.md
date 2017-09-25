Georgetown Pre-processing
================
Dominic Pearce

``` r
library(affy)
library(tidyverse)
library(Biobase)
library(frma)
library(hgu133plus2frmavecs)
library(knitr)
library(a4Base)
```

First we'll read in our .cel files in batch, normalise using frma, feature select and loess normalise - though not necessarily in that order!
---------------------------------------------------------------------------------------------------------------------------------------------

#### Create AffyBatch from .CEL files

``` r
dir_vec <- c("../data/Edinburgh\ Second\ Batch\ 170\ Samples/ff",
             "../data/Edinburgh\ Second\ Batch\ 170\ Samples/ffpe")

affybatch_lst <- lapply(dir_vec, function(dir){
                            setwd(dir)
                            affybatch <- ReadAffy()
                            setwd("../../../src")
                            affybatch
             })
```

#### Feature selection

``` r
goodcalls_lst <- lapply(affybatch_lst, function(batch){
                            #get present, marginal, absent calls
                            ap <- mas5calls(batch)
                            #for each gene check that it not called absent in less than 90% of samples
                            goodcalls <- rowSums(exprs(ap) != "A") > (ncol(ap) / 10)
                            row.names(ap)[which(goodcalls)]
             })
#retreive common probes and write out
present_common <- do.call(intersect, goodcalls_lst)
```

``` r
data.frame(present_10 = c(length(goodcalls_lst[[1]]), 
                          length(goodcalls_lst[[2]]), 
                          length(present_common)),
           material = c("FF", "FFPE", "Common")
           ) %>% knitr::kable()
```

|  present\_10| material |
|------------:|:---------|
|        41826| FF       |
|        37019| FFPE     |
|        35516| Common   |

#### 2-step normalisation - frma and loess

``` r
data(hgu133plus2frmavecs)
frma_lst <- lapply(affybatch_lst, function(batch){
                       batch_frma <- frma(batch, input.vecs = hgu133plus2frmavecs)  
                       batch_frma[present_common, ]
           })
#georgeset <- do.call(Biobase::combine, frma_lst)
#write_rds(georgeset, "~/Desktop/temp-store/georgeset-sep-frma-fselect.rds")


norm_lst <- lapply(frma_lst, function(batch){
                       mtx <- normalize.loess(exprs(batch))
                       exprs(batch) <- mtx
                       batch
           })
georgeset <- do.call(Biobase::combine, norm_lst)
#write_rds(georgeset, "~/Desktop/temp-store/georgeset-sep-frma-fselect-loess.rds")
```
