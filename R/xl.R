
#' Compute DEGs and dump with metadata, counts
#' TODO: need adjustable CPM cutoff
#' @export
deg_dump <- function(se, dir.name=".", method="QL", ...) {
    stopifnot(dir.exists(dir.name))
    A = lapply(SummarizedExperiment::assays(se), function(M) {
        cbind(data.frame(gene=rownames(M)), as.data.frame(M))
    })
    A$deg = deg.edger(se, "group", se@metadata$case, se@metadata$control, method=method, cpm.cutoff=100)
    A$gene_info = cbind(data.frame(gene=rownames(se)), SummarizedExperiment::rowData(se))
    A$sample_metadata = cbind(data.frame(name=colnames(se)), SummarizedExperiment::colData(se))
    writexl::write_xlsx(lapply(A, as.data.frame), paste0(dir.name, "/", se@metadata$comparison, ".xlsx"))
}

#' @export
deg_load <- function(xlsx) {
    comparison = gsub("[.]xlsx$", "", basename(xlsx))
    sheets = readxl::excel_sheets(xlsx)
    L = lapply(setNames(sheets, sheets), function(sheet.name) {
        as.data.frame(readxl::read_xlsx(xlsx, sheet.name))
    })
    case = unique(L$deg$case)
    control = unique(L$deg$control)
    A = lapply(L[1:(match("deg", sheets)-1)], function(af) {
        M = as.matrix(af[2:ncol(af)])
        rownames(M) = af[[1]]
        return(M)
    })
    rownames(L$deg) = L$deg$gene
    rowData=L$gene_info
    rownames(rowData) = rowData[[1]]
    rowData[[1]] = NULL
    colData=L$sample_metadata
    rownames(colData) = colData[[1]]
    colData[[1]] = NULL
    return(SummarizedExperiment::SummarizedExperiment(assays=A,
                                                      rowData=rowData,
                                                      colData=colData,
                                                      metadata=list(comparison=comparison, case=case, control=control, deg=L$deg)))
}
