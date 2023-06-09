
#' Calculate TMM assay in SummarizedExperiment
#'
#' @param se SummarizedExperiment
#' @param method Method used for normalization, e.g. TMM or RLE
#' @param log Whether to use log-CPM or CPM after normalization
#' @export
se_tmm <- function(se, method="TMM", log=FALSE) {
    dgel = edgeR::calcNormFactors(se, method)
    tmm = edgeR::cpm(dgel, log=log)
    SummarizedExperiment::assays(se)[[method]] = tmm
    return(se)
}


#' @export
se_plotPCA <- function(se, title="PCA plot", method="TMM", top=500, cpm.frac=0.25, cpm.cutoff=100, gene.selection="common", correct=FALSE, use_label=TRUE, force=0.2, text.size=4, max.overlaps=10) {
    require(ggplot2)
    dgel = edgeR::calcNormFactors(se, method)
    to_keep = rowMeans(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.frac
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    if (correct) {
        dgel = limma::removeBatchEffect(dgel, dgel$batch)
    }
    pl = limma::plotMDS(dgel, top=top, gene.selection=gene.selection, plot=FALSE)
    pf = as.data.frame(SummarizedExperiment::colData(se))
    pf$x = pl$x
    pf$y = pl$y
    if (length(unique(pf$batch))==1) {
        g = ggplot(pf, aes(x=x, y=y, label=title))
    } else {
        g = ggplot(pf, aes(x=x, y=y, color=batch, label=title)) + scale_color_brewer(palette="Set2")
    }
    g = g + geom_point()
    if (use_label) {
        g = g + ggrepel::geom_text_repel(force=force, size=text.size, max.overlaps=max.overlaps)
    }
    g = g + xlab(paste0("PC1: ", round(100 * pl$var.explained[1]), "% of variance"))
    g = g + ylab(paste0("PC2: ", round(100 * pl$var.explained[2]), "% of variance"))
    g = g + ggpubr::theme_pubr() + ggtitle(title) + theme(plot.title=element_text(size=10), legend.text=element_text(size=8))
    return(g)
}
