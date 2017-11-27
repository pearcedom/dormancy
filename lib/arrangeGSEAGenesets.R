#Given the output of GSEABase::getGmt("asdf.gmt") [GMT] file and a expressionSet [eset], organise the genesets specific to genes present in eset
arrangeGenesets <- function(GMT, eset){
    genesets <- geneIds(GMT)
    intersets <- lapply(genesets, function(set){
                            intersect(set, row.names(eset))
    })
    intersets
}
