#Two functions that take an expressionSet composed of two or more potential batches and plot expression density and mds (using 500 most variable genes) figures to determine whether to the batch differences truly exist and to what extent.
#
#eset is the expressionSet input
#batch_col is the column name of esets pData that defines the potential batch information - e.g. a column detailing whether samples were Fresh Frozen or Formaline-Fixed Parafin Embedded

print("functions are batchEffectDists() & batchEffectMDS()")

batchEffectDists <- function(eset, batch_col){
    library(Biobase)
    source("../../../../functions/mostVar.R")
    source("../../../../functions/library/mdsArrange.R")
    mtx_input <- exprs(eset)
    pheno <- pData(eset)

    #arrange
    george_mlt <- melt(mtx_input)
    george_mrg <- base::merge(george_mlt, pheno, by.x = "Var2", by.y = 0)
    george_mrg$Var2 <- factor(george_mrg$Var2, levels = unique(george_mrg$Var2[order(george_mrg[[batch_col]])]))
    #boxplot
    p_box <- ggplot(george_mrg, aes_string(x = "Var2", y = "value", fill = batch_col)) + 
        geom_boxplot() +
        theme(legend.position = 'bottom',
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    #densityplot
    p_dist <- ggplot(george_mrg, aes_string(x = "value", colour = batch_col, group = "Var2"), alpha = 0.05) + 
        geom_line(stat = 'density', alpha = 0.2) + 
        theme_pander() +    
        theme(legend.position = 'bottom',
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
        #combine
    cowplot::plot_grid(p_box, p_dist, ncol = 1)
}


batchEffectMDS <- function(eset, batch_col){
    library(Biobase)
    library(ggthemes)
    source("../../../../functions/mostVar.R")
    source("../../../../functions/library/mdsArrange.R")
    mtx_input <- exprs(eset)
    pheno <- pData(eset)
    #arrange
    mv500 <- mostVar(mtx_input, 500) %>% row.names()
    arg500 <- mdsArrange(mtx_input[mv500,]) 
    mds_input <- base::merge(arg500, pheno, by.x = 'ids', by.y = 0)
    #and plot
    ggplot(mds_input, aes_string("x", "y", colour = batch_col)) + 
        geom_point() + 
        theme_pander() + 
        theme(legend.position = 'bottom')
}
