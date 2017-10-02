Clincal summaries
================
Dominic Pearce

``` r
library(readr)
library(knitr)
library(Biobase)
```

Georgetown dormancy data summaries
----------------------------------

``` r
georgeset <- read_rds("../output/final-georgeset-sep-frma-fselect-loess-clin-cb.Rds")

columns_g <- c("is_dormant", "timepoint", "Age", "Surgery.Performed", "T", "N", "M", "ER", "Diagnostic.CB.Grade", "HER2.IHC", "HER2.FISH", "Overall.HER2", "No.of.Pos.Nodes", "Chemo", "Radiotherapy", "ET.Drugs", "Any.Recurrence", "Alive.Dead", "is_freshfrozen")

for(x in columns_g){
                   dfr <- data.frame(table(pData(georgeset)[,x]))
                   colnames(dfr) <- c(x, "Freq")
                   kable(dfr)
                } 
```

Edinburgh dormancy data summaries
---------------------------------

``` r
edset <- read_rds("../../edinburgh/output/dorm-v4.rds")

columns_e <- c("dorm.group_v4", "time.point_3cat", "age", "T", "N", "M", "ER", "Grade", "HER2", "HER2.FISH", "drug.type", "proteomics", "recur.status", "Death")

for(x in columns_e){
           dfr <- data.frame(table(pData(edset)[,x]))
           colnames(dfr) <- c(x, "Freq")
           print(kable(dfr))
}
```

    ## 
    ## 
    ## dorm.group_v4    Freq
    ## --------------  -----
    ## D                 116
    ## ND                 61
    ## 
    ## 
    ## time.point_3cat    Freq
    ## ----------------  -----
    ## 1                    60
    ## 2                    52
    ## 4                    65
    ## 
    ## 
    ## age        Freq
    ## --------  -----
    ##             116
    ## 53            1
    ## 56            2
    ## 57            3
    ## 60            1
    ## 61            1
    ## 63            2
    ## 64            1
    ## 65            2
    ## 67            1
    ## 68            2
    ## 69            1
    ## 71            1
    ## 72            3
    ## 73            3
    ## 74            2
    ## 75            6
    ## 76            1
    ## 78            5
    ## 79            4
    ## 80            1
    ## 81            1
    ## 83            1
    ## 84            2
    ## 85            2
    ## 87            3
    ## 89            1
    ## Unknown       8
    ## 
    ## 
    ## T          Freq
    ## --------  -----
    ## T1           24
    ## T2           78
    ## T3            9
    ## T4           10
    ## T4b          21
    ## T4B           8
    ## T4D           4
    ## Unknown      20
    ## 
    ## 
    ## N          Freq
    ## --------  -----
    ## N0           98
    ## N1           47
    ## N3            1
    ## NO            2
    ## NX            3
    ## Unknown      23
    ## 
    ## 
    ## M          Freq
    ## --------  -----
    ## M0          139
    ## M1            7
    ## MX            5
    ## Unknown      23
    ## 
    ## 
    ## ER    Freq
    ## ---  -----
    ## 6        3
    ## 7       29
    ## 8      145
    ## 
    ## 
    ## Grade      Freq
    ## --------  -----
    ## 1            19
    ## 2           101
    ## 3            41
    ## Unknown      16
    ## 
    ## 
    ## HER2                Freq
    ## -----------------  -----
    ## 0                     14
    ## 1+                    12
    ## 2+                    37
    ## 3+                    22
    ## neg                   31
    ## Neg (Gene Level)       4
    ## negative              24
    ## Negative               3
    ## positive               2
    ## unknown                2
    ## Unknown                7
    ## Unknown               13
    ## Unknown(neg)           3
    ## 
    ## 
    ## HER2.FISH    Freq
    ## ----------  -----
    ##                36
    ## Negative       55
    ## Positive       22
    ## unknown        22
    ## Unknown        39
    ## 
    ## 
    ## drug.type      Freq
    ## ------------  -----
    ##                 151
    ## Anastrazole      20
    ## Tamoxifen         6
    ## 
    ## 
    ## proteomics    Freq
    ## -----------  -----
    ## 0              130
    ## 1               47
    ## 
    ## 
    ## recur.status    Freq
    ## -------------  -----
    ## 0                120
    ## 1                 42
    ## Unknown           15
    ## 
    ## 
    ## Death      Freq
    ## --------  -----
    ## 0            86
    ## 1            76
    ## Unknown      15
