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

lapply(columns_g, function(x){
                   dfr <- data.frame(table(pData(georgeset)[,x]))
                   colnames(dfr) <- c(x, "Freq")
                   kable(dfr)
                }) 
```

    ## [[1]]
    ## 
    ## 
    ## is_dormant    Freq
    ## -----------  -----
    ## FALSE           32
    ## TRUE            71
    ## 
    ## [[2]]
    ## 
    ## 
    ## timepoint       Freq
    ## -------------  -----
    ## diagnosis         51
    ## long-term         50
    ## on-treatment      36
    ## 
    ## [[3]]
    ## 
    ## 
    ## Age    Freq
    ## ----  -----
    ## 57        2
    ## 59        1
    ## 60        3
    ## 61        4
    ## 63        7
    ## 64        4
    ## 65        2
    ## 67        2
    ## 69        7
    ## 71        3
    ## 72        7
    ## 73       18
    ## 75       14
    ## 76        5
    ## 78       13
    ## 79        3
    ## 80        1
    ## 81        2
    ## 84       20
    ## 85        5
    ## 88        4
    ## 89        5
    ## 92        2
    ## 
    ## [[4]]
    ## 
    ## 
    ## Surgery.Performed    Freq
    ## ------------------  -----
    ## mast                    9
    ## Mast                    1
    ## no surgery              9
    ## WLE                   113
    ## 
    ## [[5]]
    ## 
    ## 
    ## T      Freq
    ## ----  -----
    ## T1       17
    ## T2       48
    ## T3       17
    ## T4       22
    ## T4b      10
    ## T4d       4
    ## 
    ## [[6]]
    ## 
    ## 
    ## N     Freq
    ## ---  -----
    ## N0      77
    ## N1      36
    ## N2      12
    ## 
    ## [[7]]
    ## 
    ## 
    ## M     Freq
    ## ---  -----
    ## M0     116
    ## M1       2
    ## 
    ## [[8]]
    ## 
    ## 
    ## ER    Freq
    ## ---  -----
    ## 7       14
    ## 8      118
    ## 
    ## [[9]]
    ## 
    ## 
    ## Diagnostic.CB.Grade    Freq
    ## --------------------  -----
    ## 1                         7
    ## 2                        86
    ## 3                        41
    ## 
    ## [[10]]
    ## 
    ## 
    ## HER2.IHC    Freq
    ## ---------  -----
    ## 0             30
    ## 1+            23
    ## 2+            39
    ## 3+             7
    ## 
    ## [[11]]
    ## 
    ## 
    ## HER2.FISH    Freq
    ## ----------  -----
    ## neg            34
    ## pos             5
    ## 
    ## [[12]]
    ## 
    ## 
    ## Overall.HER2    Freq
    ## -------------  -----
    ## neg               87
    ## pos               12
    ## 
    ## [[13]]
    ## 
    ## 
    ## No.of.Pos.Nodes    Freq
    ## ----------------  -----
    ## 0                    62
    ## 1                    17
    ## 2                    11
    ## 3                     6
    ## 4                     3
    ## 6                     4
    ## 8                    15
    ## 9                     4
    ## 10                    5
    ## 
    ## [[14]]
    ## 
    ## 
    ## Chemo    Freq
    ## ------  -----
    ## No         80
    ## Yes        17
    ## 
    ## [[15]]
    ## 
    ## 
    ## Radiotherapy    Freq
    ## -------------  -----
    ## No                21
    ## Yes               76
    ## 
    ## [[16]]
    ## 
    ## 
    ## ET.Drugs         Freq
    ## --------------  -----
    ## Let               114
    ## Let-Anas            8
    ## Let-Anas-Exem       1
    ## Let-Exem            3
    ## Let-Tam             3
    ## Let-Tam-Let         5
    ## 
    ## [[17]]
    ## 
    ## 
    ## Any.Recurrence                                 Freq
    ## --------------------------------------------  -----
    ## N/A                                               7
    ## No                                               93
    ## Yes                                              20
    ## Yes (Liver mets at diagnosis and recurrent)       2
    ## 
    ## [[18]]
    ## 
    ## 
    ## Alive.Dead    Freq
    ## -----------  -----
    ## Alive           44
    ## Dead            68
    ## 
    ## [[19]]
    ## 
    ## 
    ## is_freshfrozen    Freq
    ## ---------------  -----
    ## FALSE              102
    ## TRUE                68

Edinburgh dormancy data summaries
---------------------------------

``` r
edset <- read_rds("../../edinburgh/output/dorm-v4.rds")

colnames(pData(edset))
```

    ##  [1] "study.type"          "dorm.group"          "dorm.group_v2"      
    ##  [4] "dorm.group_v3"       "dg_v3_evidence"      "dorm.group_v4"      
    ##  [7] "dg_v4_criteria"      "v4_uss"              "v4_pmarker"         
    ## [10] "PAM50"               "timetolastUSS"       "USS_class"          
    ## [13] "prog_class"          "prog_status"         "prog_status_v3"     
    ## [16] "time_to_prog"        "patientID"           "patient.ID"         
    ## [19] "patient.no"          "age"                 "biopsy.no"          
    ## [22] "sampling.type"       "days_newinfo"        "resp.group"         
    ## [25] "recur.status"        "recur.type"          "timetorec"          
    ## [28] "Death"               "BCSpecific"          "Osurvival"          
    ## [31] "DFS"                 "T"                   "N"                  
    ## [34] "M"                   "Grade"               "ER"                 
    ## [37] "HER2"                "HER2.FISH"           "HER_myfinaldecision"
    ## [40] "onTMA"               "drug.change"         "drug.type"          
    ## [43] "proteomics"          "ID_D120days"         "time.point"         
    ## [46] "ID_D120days_3cat"    "time.point_3cat"

``` r
columns_e <- c("dorm.group_v4", "time.point_3cat", "age", "T", "N", "M", "ER", "Grade", "HER2", "HER2.FISH", "drug.type", "proteomics", "recur.status", "Death")

obs <- lapply(columns_e, function(x){
           dfr <- data.frame(table(pData(edset)[,x]))
           colnames(dfr) <- c(x, "Freq")
           print(kable(dfr))
                })
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
