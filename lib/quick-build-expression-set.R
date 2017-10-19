#Build expressionSet from matched expression and clinical matrices
buildeSet <- function(xpr, pheno){
    library(Biobase)
    library(testthat)
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
