
.proc.go.result <- function(obj) {
    if (is.null(obj)) {
        return(data.frame())
    }
    if ("result" %in% slotNames(obj)) {
        df = slot(obj, "result")
        if (nrow(df) == 0) {
            return(data.frame())
        } else {
            return(df)
        }
    } else {
        return(data.frame())
    }
}
#' Run clusterProfiler::enrichGO with better options and defaults
#'
#' @param gene List of genes
#' @param OrgDb OrgDb e.g. org.Hs.eg.db::org.Hs.eg.db or org.Mm.eg.db::org.Mm.eg.db
#' @param keyType Set default keyType to SYMBOL instead of ENTREZID for clusterProfiler::enrichGO
#' @param ont Set default ontology to ALL for clusterProfiler::enrichGO
#' @param minGSSize Minimum gene set size for a GO term
#' @param maxGSSize Maximum gene set size for a GO term
#' @export
enrichGO <- function(gene, OrgDb, keyType="SYMBOL", ont="ALL", minGSSize=10, maxGSSize=100, ...) {
    ego = clusterProfiler::enrichGO(gene=gene, OrgDb=OrgDb, keyType=keyType, ont=ont, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
    return(.proc.go.result(ego))
}


#' @export
top.go <- function(se, fdr.cutoff=0.05, log2FC.cutoff=0, cpm.cutoff=7.5, ...) {
    rd = SummarizedExperiment::rowData(se)
    cd = SummarizedExperiment::colData(se)
    md = S4Vectors::metadata(se)
    deg = md$deg[(md$deg$FDR < fdr.cutoff) & (md$deg$log2FC > log2FC.cutoff) & (md$deg$logCPM >= cpm.cutoff),]
    bg = unique(md$deg$gene)
    enrichGO(gene=deg$gene, OrgDb=org.Hs.eg.db::org.Hs.eg.db, universe=bg, ...)
}
