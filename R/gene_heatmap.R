

#' Plot a heatmap for the following genes.
#' SE must be sorted by group.
#' @export
gene_heatmap <- function(se, genes=NULL, max.ngene=50, assay="TMM", group="group", column_title="Heatmap of genes", ux=4, balance.deg=FALSE) {
    ### get data
    rd = SummarizedExperiment::rowData(se)
    cd = SummarizedExperiment::colData(se)
    md = S4Vectors::metadata(se)
    M = SummarizedExperiment::assays(se)[[assay]]
    if (!is.null(genes)) {
        M = M[intersect(genes, rownames(M)),]
    }
    deg = md$deg[md$deg$gene %in% rownames(M),] ## keep order
    m.ord = rep(Inf, nrow(M))
    if (balance.deg) {
        deg.pos = deg[deg$log2FC > 0,]
        deg.neg = deg[deg$log2FC < 0,]
### reorder M based on deg literal order
        m.ord[match(deg.pos$gene, rownames(M))] = 1:nrow(deg.pos)
        m.ord[match(deg.neg$gene, rownames(M))] = 1:nrow(deg.neg)
    } else {
        m.ord[match(deg$gene, rownames(M))] = 1:nrow(deg)
    }
    M = M[order(m.ord),]
    if (nrow(M) > max.ngene) {
        M = M[1:max.ngene,]
    }
    ### after DEG subset, order by log2FC instead of average expression
    M = M[order(deg[rownames(M),]$log2FC),]
    ### pseudobulk
    P = make_pseudobulk(cd[[group]])
    MP = vclip(M %*% P %*% diag(1/colSums(P)), min=0)
    colnames(MP) = colnames(P)
    ### transpose
    M = scale(t(M))
    MP = t(MP)
    ### set pvalues
    pvalues = setNames(rep(1, ncol(M)), colnames(M))
    pvalues[deg$gene] = deg$FDR
    ### order based on pseudobulk in #1
    M = M[,order(-MP[1,colnames(M)])]
    MP = MP[,colnames(M)]
    pvalues = pvalues[colnames(M)]
    ### annotation
    ta = ComplexHeatmap::HeatmapAnnotation(sig=ComplexHeatmap::anno_text(convert_pvalue_to_star(pvalues), just="center",
                                                                         rot=0, location=grid::unit(0.5, 'npc'),
                                                                         gp = grid::gpar(fontsize=5)
                                                                         ),
                                           show_legend=FALSE,
                                           show_annotation_name=TRUE,
                                           annotation_name_gp= grid::gpar(fontsize = 3.5),
                                           annotation_name_side="right",
                                           `Mean Expression`=ComplexHeatmap::anno_barplot(log1p(t(MP))/log(10),
                                                                             beside=TRUE, attach=TRUE, #baseline="min",
                                                                                        #gp = grid::gpar(fill="#CCCCCC"),
                                                                                                axis_param=list(gp=grid::gpar(fontsize=3), side="left"),
                                                                             border=FALSE, height=grid::unit(0.5, "cm")))
    ra = ComplexHeatmap::rowAnnotation(`Counts per sample`=ComplexHeatmap::anno_barplot(colSums(SummarizedExperiment::assays(se)[["counts"]]), border=FALSE, height=grid::unit(0.75, "cm"), axis_param=list(gp=grid::gpar(fontsize=3))),
                                       annotation_name_gp= grid::gpar(fontsize = 3.5),
                                       show_legend=FALSE,
                                       show_annotation_name=TRUE)

### color
    col = circlize::colorRamp2(c(-1.5, #min(-1.5, quantile(M, 0.1, na.rm=TRUE)),
                                 0,
                                 1.5), #max(1.5, quantile(M, 0.9, na.rm=TRUE))),
                               c("blue", "white", "red"))
    H = autoHeatmap(M, ux=ux, dimname_fontsize=6,
                    top_annotation=ta,
                    right_annotation=ra,
    #                col=col,
                    row_split=cd$group,
                    row_title_gp=grid::gpar(fontsize=8, font=2),
                    column_title_gp=grid::gpar(fontsize=8, font=2),
                    cluster_row_slices=FALSE, sort=NULL,
                    heatmap_legend_param=list(title="Scaled gene expression",
                                              title_gp = grid::gpar(fontsize = 3.5, font = 2),
                                              legend_gp = grid::gpar(fontsize = 5),
                                              labels_gp = grid::gpar(fontsize = 5),
                                              grid_height = grid::unit(4, "mm"), grid_width = grid::unit(3, "mm")
                                             ),
                    column_title=column_title,
                    row_title_rot = 0,
                    show_row_names=FALSE,
                    #row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                    #column_title_gp = grid::gpar(fontsize = 3, fontface = "bold")
                   )
    w = (ncol(M) + 18) * ux * 0.03937
    h = (nrow(M) + 8) * ux * 0.03937
    attr(H, "mib_w") = w
    attr(H, "mib_h") = h
    return(H)
}
