
#' Add a fancy barplot for a specific gene
#' @export
gene_barplot <- function(se, gene, assay="TMM", group="group", add=c("jitter"), method="t.test") {
    cd = SummarizedExperiment::colData(se)
    cd = as.data.frame(cd)
    M = SummarizedExperiment::assays(se)[[assay]]
    cd[[gene]] = M[gene,]
    ggpubr::ggboxplot(cd, x=group, y=gene, add=add) + ggpubr::stat_compare_means(method=method)
}
