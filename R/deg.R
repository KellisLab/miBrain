## DEG methods from https://github.com/KellisLab/benj

#' Run edgeR on a SummarizedExperiment
#'
#' @param se SummarizedExperiment to run
#' @param pathology Column in colData(se) that contains case and control
#' @param case Case string within pathology column. Example "AD"
#' @param control Control string within pathology column. Example "control"
#' @param covariates Vector of covariates to add to design matrix
#' @param method Method used for edgeR. Either LRT or QL
#' @param assay Assay from se to utilize
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.frac Fraction of samples passing CPM cutoff for filter
#' @export
deg.edger <- function(se, pathology, case, control, covariates=c(),
                      method="LRT", assay=NULL, cpm.cutoff=10, cpm.frac=0.25,
                      filter_only_case_control=FALSE) {
    if (filter_only_case_control) {
        se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
    }

    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)$counts
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    if (S4Vectors::ncol(X) > 1000) {
        warning("More than 1000 samples are being called to a function under the assumption data passed in is pseudobulk")
    }

    cd = SummarizedExperiment::colData(se)
    if (!is.character(cd[[pathology]]) | !is.factor(cd[[pathology]])) {
        warning(paste0("Converting ", pathology, " to factor"))
    }
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    to_keep = Matrix::rowMeans(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.frac
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    dgel = edgeR::calcNormFactors(dgel, method="TMM")
    design = model.matrix(as.formula(paste0("~0 + ", paste0(c(pathology, covariates), collapse=" + "))),
                          data=cd)
    colnames(design) = make.names(colnames(design)) ### need spaces to be removed
    contrasts = limma::makeContrasts(contrasts=paste0(make.names(paste0(pathology, case)),
                                                      "-",
                                                      make.names(paste0(pathology, control))),
                                     levels=design)
    if (!(method %in% c("QL", "LRT"))) {
        method = "QL"
        warning(paste0("Setting method to ", method))
    }
    if (method == "QL") {
        dgel = edgeR::estimateDisp(dgel, design, robust=TRUE)
        fit <- edgeR::glmQLFit(dgel, design, robust=TRUE)
        res <- edgeR::glmQLFTest(fit, contrast=contrasts)
    } else if (method == "LRT") {
        dgel = edgeR::estimateGLMCommonDisp(dgel, design)
        dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
        fit = edgeR::glmFit(dgel, design, robust=TRUE)
        res = edgeR::glmLRT(fit, contrast=contrasts)
    }
    tbl = edgeR::topTags(res, n=Inf, sort.by="none")$table
    colnames(tbl) = gsub("^logFC$", "log2FC", colnames(tbl))
    tbl$gene = rownames(tbl)
    tbl$method = method
    tbl$case = case
    tbl$control = control
    tbl$logCPM = edgeR::cpm(Matrix::rowSums(X), log=TRUE)[tbl$gene,]
    return(tbl[order(tbl$FDR),])
}


#' Run DESeq2 on pseudobulked data
#'
#' @param se SummarizedExperiment data, expected pseudobulk
#' @param pathology Column in colData(se) from which to study differential changes
#' @param case Case value within pathology column
#' @param control Control value within pathology column
#' @param covariates Covariates to pass to DESeq2
#' @param assay Assay to perform differential changes upon. Default is "counts"
#' @param filter_only_case_control logical telling whether to filter sce to only case and control values within pathology
#' @export
deg.deseq2 <- function(se,
                       pathology,
                       case,
                       control,
                       covariates=c(),
                       assay=NULL,
                       filter_only_case_control=FALSE) {
    if (filter_only_case_control) {
        se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
    }
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)$counts
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    if (S4Vectors::ncol(X) > 1000) {
        warning("More than 1000 samples are being called to a function under the assumption data passed in is pseudobulk")
    }
    formula = paste0("~", c(pathology, covariates), collapse=" + ")
    dds = DESeq2::DESeqDataSetFromMatrix(
                      X,
                      colData=SummarizedExperiment::colData(se),
                      design=as.formula(formula))
    out = DESeq2::DESeq(dds)
    df = DESeq2::results(out, contrast=c(pathology, case, control))
    df = df[!is.na(df$padj),]
    df = df[order(df$padj),]
    df$gene = rownames(df)
    colnames(df) = msub(c("^log2FoldChange$", "^padj$"),
                        c("log2FC", "FDR"), colnames(df))
    df$case = case
    df$control = control
    df$logCPM = edgeR::cpm(Matrix::rowSums(X), log=TRUE)[df$gene,]
    return(as.data.frame(df[order(df$FDR),]))
}
