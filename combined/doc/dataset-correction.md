Dataset correction
================
Dominic Pearce

#### //TODO : I've ignored recurrence data below but this needs to be integrated into the integration as well

#### //TODO : When using the updated edset data, there'll likely have to be a couple of bits and pieces to fix from Cigdem's output

#### //TODO use better idReplace version (probe deicison by highest summed total)

#### //TODO recheck dormancy assignment consistency between datasets once xpr-based assigments have been made for georgeset

#### //TODO The final reassignment step should really be performed in clinical-curation.rmd

``` r
library(tidyverse)
library(Biobase)
library(sva)
library(testthat)
source("/Volumes/igmm/sims-lab/Dominic/functions/idReplace.R")
source("/Volumes/igmm/sims-lab/Dominic/functions/convert-id-to-gene-symbol-with-biomart.R")
```

#### Now we have the two dormancy datasets - Edinburgh & Georgetown - pre-processed and ready to analyse. In addition to their independent analysis, we'll also want to pool all samples into a single dataset for maximum power. **However**, having been pre-processed months apart, in different countries and by different people, before combining we first need to make sure their clinical data is analogous.

#### Ideally we want to keep as much information as possible but ultimately there will be fields in one eset that aren't found in the other.

``` r
ggset <- read_rds("../../georgetown/output/final-georgeset-sep-frma-fselect-loess-clin-cb.Rds")
edset <- read_rds("../../edinburgh/output/dorm-v4.rds")
edset$platform <- "illumina_humanht_12"
edset$xpr_id <- colnames(edset)

ggsub <- pData(ggset)[, c("sample_id", "patient", "is_dormant", "timepoint", "days_treated", "biopsy", "Age", "T", "N", "M", "ER", "Overall.HER2", "platform", "xpr_id")]
edsub <- pData(edset)[, c("patientID", "patient.no", "dorm.group_v4", "time.point_3cat", "days_newinfo", "biopsy.no", "age", "T", "N", "M", "ER", "HER_myfinaldecision", "platform", "xpr_id")]
new_cols <- c("sample_id", "patient_id", "is_dormant", "timepoint", "days_treated", "biopsy_no", "age", "T", "N", "M", "ER_allred", "is_her2", "platform", "xpr_id")
colnames(ggsub) <- new_cols
colnames(edsub) <- new_cols

edsub$is_dormant <- edsub$is_dormant == "D"
edsub$timepoint <- ifelse(edsub$timepoint == "1", "diagnosis", 
                          ifelse(edsub$timepoint == "2", "on-treatment", "long-term"))
edsub$is_her2 <- ifelse(edsub$is_her2 == "1", TRUE, 
                        ifelse(edsub$is_her2 == "0", FALSE, NA))

ggsub$is_her2 <- ifelse(ggsub$is_her2 == "pos", TRUE, 
                        ifelse(ggsub$is_her2 == "neg", FALSE, NA))
```

#### Convert georgtown xpr to official gene symbols

``` r
gghgnc <- idReplace(exprs(ggset), "affy_hg_u133_plus_2", "hgnc_symbol")
```

#### Now with our datasets made comparable we combine phenotypic and assay data as a single expression set `dormset`

``` r
#Combine expression
dorm_xpr <- merge(exprs(edset), gghgnc, by = 0)
row.names(dorm_xpr) <- dorm_xpr$Row.names
dorm_xpr$Row.names <- NULL

#Combine clinical data
dorm_pheno <- rbind(edsub, ggsub)

#Determine shared samples. For this we first strip any letters that may be in patient ids
dorm_pheno$patient <- gsub("[A-Z]", "", dorm_pheno$patient)
dorm_pheno$sample_id <- gsub("[A-Za-z]", "", dorm_pheno$sample_id)

shared_samples <- names(which(table(dorm_pheno$sample_id) != 1))
shared_patients <- unique(gsub("\\-[0-9]", "", shared_samples))
dorm_pheno$is_dataset_replicate <- dorm_pheno$sample_id %in% shared_samples

test_that("dimnames are correct between xpr and pheno", {
              expect_identical(row.names(dorm_pheno), colnames(dorm_xpr))
                        })
```

#### So the shared sample n and additional sample n are... (pending further clinical information)

``` r
length(shared_samples) #57 shared samples
```

    ## [1] 68

``` r
length(shared_patients) #24 shared patients
```

    ## [1] 24

``` r
addset <- ggset[, which(!ggset$sample_id %in% shared_samples)]
sharedset <- ggset[, which(ggset$sample_id %in% shared_samples)]
#additional dormant/desensitised samples
table(addset$is_dormant, exclude = NULL)
```

    ## 
    ## FALSE  TRUE  <NA> 
    ##     6    37    14

``` r
addset$patient %>% unique() %>% length()
```

    ## [1] 24

``` r
#additional dormant patients
addset[, addset$is_dormant]$patient %>% unique() %>% length()
```

    ## [1] 16

``` r
#additional desensitised patients
addset[, !addset$is_dormant]$patient %>% unique() %>% length()
```

    ## [1] 4

``` r
#Build expressionSet
buildeSet <- function(xpr, pheno){
    test_that("dimnames are correct between xpr and pheno", {
                  expect_identical(row.names(pheno), colnames(xpr))
                        })
    metadata <- data.frame(labelDescription = colnames(pheno), row.names = colnames(pheno))
    phenoData <- new("AnnotatedDataFrame", data = pheno, varMetadata = metadata)
    ExpressionSet(
                  assayData = as.matrix(xpr),
                  phenoData = phenoData
                  )
}
dormset <- buildeSet(dorm_xpr, dorm_pheno)
```

#### Now we check for batch differences, correct using ComBat and validate our approach using the samples which are replicated between the two datasets.

#### *PRE-ComBat*

``` r
source("../../lib/plot-batch-effects.R")
```

    ## [1] "functions are batchEffectDists() & batchEffectMDS()"

``` r
batchEffectDists(dormset, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(dormset, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

#### *POST-ComBat*

``` r
dorm_cb <- ComBat(exprs(dormset), batch = dormset$platform)
```

    ## Standardizing Data across genes

``` r
dormset_cb <- buildeSet(dorm_cb, pData(dormset))

batchEffectDists(dormset_cb, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(dormset_cb, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-2.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(dormset_cb, "timepoint")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-3.png" style="display: block; margin: auto;" />

### Replicates & replicate correlations

#### Replicate-only MDS

``` r
batchEffectMDS(dormset_cb[, dormset_cb$is_dataset_replicate], "platform") + 
    geom_point(colour = "WHITE") + 
    geom_text(aes(label = sample_id))
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

#### ComBat appears to have worked well, with replicate samples matching closely. The final step is to remove the unwanted replicate samples and re-batch correct using only the samples that will part of the final analysis

#### Because all Edinburgh samples were fresh frozne, any Georgetown replicates that were FFPE can simply be removed.

``` r
ffpe_to_remove <- sharedset$sample_id[!sharedset$is_freshfrozen]
```

#### For the FF samples however, we'll calculate their intra-dataset correlations, average this across a patient and select the samples from the dataset that has the highest patient-specific correlation for each individual set of patient samples. Or more simply, we won't just keep samples with the best correlations, we'll select patients who's samples average the best correlation.

#### This will therefore avoid selecting samples 123-1 & 123-2 from the Edinburgh dataset and 123-3 from the Georgetown dataset, and hence avoid correcting within a patient sample set.

``` r
ff_dupes <- sharedset$sample_id[sharedset$is_freshfrozen]

dupeCors <- function(platform){
    eset <- dormset[, dormset$platform == platform]
    xpr <- exprs(eset)
    colnames(xpr) <- eset$sample_id
    set_cor <- cor(xpr)
    diag(set_cor) <- NA
    set_cormeans <- colMeans(set_cor, na.rm = TRUE)
    set_cormeans[names(set_cormeans) %in% ff_dupes]
}

sharedcor_dfr <- data.frame(georgetown = dupeCors("affy_hg_u133_plus_2"), 
                            edinburgh = dupeCors("illumina_humanht_12"), 
                            patient = gsub("\\-[0-9]", "", ff_dupes))

shared_decision <- sapply(unique(sharedcor_dfr$patient), function(x){
                              dfr <- sharedcor_dfr[sharedcor_dfr$patient == x, 1:2]
                              sel <- which.max(colMeans(dfr))
                              c("georgetown", "edinburgh")[sel]
                            })
names(shared_decision) <- unique(sharedcor_dfr$patient)
shared_decision %>% table()
```

    ## .
    ##  edinburgh georgetown 
    ##          4         15

#### It's evident that the &gt;75% of the patient sets will be taken from the georgetown data. Therefore, because the Georgetown data has already been reduced by removing FFPE samples, we'll keep all patient sets from Georgetown, to avoid as best we can unbalancing our batches during ComBat.

``` r
dormset_unique <- dormset[,-c(
                             which(dormset$platform == "affy_hg_u133_plus_2" &
                                   dormset$sample_id %in% ffpe_to_remove),
                             which(dormset$platform == "illumina_humanht_12" &
                                   dormset$sample_id %in% ff_dupes)
                                   )]

dorm_cb_unique <- ComBat(exprs(dormset_unique), batch = dormset_unique$platform)
```

    ## Standardizing Data across genes

``` r
dormset_cb_unique <- buildeSet(dorm_cb_unique, pData(dormset_unique))

batchEffectDists(dormset_cb_unique, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(dormset_cb_unique, "platform")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-14-2.png" style="display: block; margin: auto;" />

``` r
batchEffectMDS(dormset_cb_unique, "timepoint")
```

<img src="dataset-correction_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-14-3.png" style="display: block; margin: auto;" />

#### Fix any weirdness that may have crept into the clinical info, order the in a sensible way and then write

#### Misc tidying up steps

``` r
dormset_out1 <- dormset_cb_unique[, order(dormset_cb_unique$sample_id,
                                         dormset_cb_unique$days_treated)]

dormset_out2 <- dormset_out1[, !is.na(dormset_out1$is_dormant)]

patient_vec <- unique(dormset_out2$patient)

#despite not having days_treated info for all samples, if the last sample that does have days_treated data is long-term, then any subsetquent biopsies must also be long-term
dormset_out3 <- dormset_out2
pData(dormset_out3) <- lapply(patient_vec, function(patient){
                           dfr <- pData(dormset_out2)[dormset_out2$patient == patient, ] 
                           lts <- which(dfr$timepoint == 'long-term')
                           nas <- which(is.na(dfr$timepoint))
                           if(length(nas) >= 1 & 
                              length(lts) >= 1 & 
                              isTRUE(nas[1] > lts[length(lts)])){
                               dfr[nas, 'timepoint'] <- 'long-term'
                           }
                           dfr
}) %>% do.call(rbind, .)


test_that("we've only changed those specific values in $timepoint", {
              expect_identical(dim(dormset_out2), dim(dormset_out3))
              expect_identical(pData(dormset_out2)[, 1:3], 
                               pData(dormset_out3)[, 1:3])
})

#designate early/late biopsies if two fall within the same timepoint - i.e. if 
#there are two long-term samples from the same patient
dormset_out3$timepoint_split <- lapply(patient_vec, function(patient){
           dfr <- pData(dormset_out3)[dormset_out3$patient == patient, ] 
           tp_tbl <- table(dfr$timepoint)
           multi_vec <- lapply(unique(dfr$timepoint), function(tp){
                                   if(is.na(tp)){
                                       NA
                                   }
                                   else if(tp_tbl[tp] > 1){
                                       if(tp == 'on-treatment'){
                                        LETTERS[1:tp_tbl[tp]]
                                       } else if(tp == 'long-term')
                                        LETTERS[tp_tbl[tp]:1]
                                   } else {
                                       NA
                                   }
                          }) %>% unlist()
           multi_vec
}) %>% do.call(c, .)

#check against set2 again
test_that("we've only changed those specific values in $timepoint", {
              expect_identical(dim(dormset_out2), dim(dormset_out3))
              expect_identical(pData(dormset_out2)[, 1:3], 
                               pData(dormset_out3)[, 1:3])
})

#patients 298 (D>AR), 341 (AR>D), 400 (AR>D) and 347 (D>AR) need to be reclassified according to xpr measurements
dormset_out3$is_dormant <- sapply(dormset_out3$patient, function(x)
       ifelse(x == 298, FALSE, dormset_out3$is_dormant[dormset_out3$patient == x])
       )
dormset_out3$is_dormant <- sapply(dormset_out3$patient, function(x)
       ifelse(x == 341, TRUE, dormset_out3$is_dormant[dormset_out3$patient == x])
       )
dormset_out3$is_dormant <- sapply(dormset_out3$patient, function(x)
       ifelse(x == 400, TRUE, dormset_out3$is_dormant[dormset_out3$patient == x])
       )
dormset_out3$is_dormant <- sapply(dormset_out3$patient, function(x)
       ifelse(x == 347, FALSE, dormset_out3$is_dormant[dormset_out3$patient == x])
       )
```

#### Write

``` r
#write_rds(dormset_out3, "../output/dormset.rds")
```
