
#' @export
PlotPCA <- function(se, title=NULL, group="group", method="TMM", top=500, cpm.frac=0.25, cpm.cutoff=100, gene.selection="common", correct=FALSE, use_label=FALSE, force=0.2, text.size=4, max.overlaps=10) {
    require(ggplot2)
    title = ifelse(is.null(title), S4Vectors::metadata(se)$comparison, title)
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
    if (is.null(group)) {
        g = ggplot(pf, aes(x=x, y=y, label=title))
    } else {
        g = ggplot(pf, aes(x=x, y=y, color=group, label=title)) + scale_color_brewer(palette="Set2")
    }
    g = g + geom_point()
    if (use_label) {
        g = g + ggrepel::geom_text_repel(force=force, size=text.size, max.overlaps=max.overlaps)
    }
    g = g + xlab(paste0("PC1: ", round(100 * pl$var.explained[1]), "% of variance"))
    g = g + ylab(paste0("PC2: ", round(100 * pl$var.explained[2]), "% of variance"))
    g = g + ggpubr::theme_pubr() + ggtitle(title) + theme(plot.title=element_text(size=11), legend.text=element_text(size=4.5), legend.title=element_text(size=8))
    return(g)
}
