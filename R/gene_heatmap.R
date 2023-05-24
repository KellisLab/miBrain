
#' @export
gene_heatmap <- function(se, genes=NULL, assay="TMM", scale=TRUE, log=TRUE, ...) {
    rd = SummarizedExperiment::rowData(se)
    cd = SummarizedExperiment::colData(se)
    md = S4Vectors::metadata(se)
    M = SummarizedExperiment::assays(se)[[assay]]
    if (log) {
        stopifnot(min(M) >= 0)
        print("Logging")
        M = base::log1p(M)
    }
    if (scale) {
        print("Scaling")
        M = base::scale(M)
    }
    if (!is.null(genes)) {
        M = M[genes,]
        rd = rd[genes,]
    }
    autoHeatmap(M, sort=NULL, ...)
}
