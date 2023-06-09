
#' Set default heatmap options, modified from Carles Boix
#'
#' @param htsc Padding ratio
#' @export
set_ht_opt <- function(htsc=2.5) {
    ComplexHeatmap::ht_opt(
        heatmap_row_names_gp = grid::gpar(fontsize = 5),
        heatmap_column_names_gp = grid::gpar(fontsize = 5),
        heatmap_row_title_gp = grid::gpar(fontsize = 5.5, fontface="bold"),
        heatmap_column_title_gp = grid::gpar(fontsize = 5.5, fontface="bold"),
        legend_title_gp = grid::gpar(fontsize = 5.5, font=2),
        legend_labels_gp = grid::gpar(fontsize = 5),
        legend_grid_height = grid::unit(2.5, 'mm'),
        legend_grid_width = grid::unit(2.5, 'mm'),
        DENDROGRAM_PADDING = grid::unit(.5 / htsc, 'mm'),
        DIMNAME_PADDING = grid::unit(1 / htsc, 'mm'),
        COLUMN_ANNO_PADDING = grid::unit(1 / htsc, 'mm'),
        ROW_ANNO_PADDING = grid::unit(1 / htsc, 'mm'),
        HEATMAP_LEGEND_PADDING = grid::unit(2 / htsc, 'mm'),
        ANNOTATION_LEGEND_PADDING = grid::unit(2 / htsc, 'mm'))
}

savePlot <- function(p, pltprefix, w=7, h=7, dpi=600, ...){
  p <- substitute(p)
  pdf(paste0(pltprefix, ".pdf"), width=w, height=h)
  eval(p)
  dev.off()

  png(paste0(pltprefix, ".png"), res=dpi, width=w, height=h, units="in")
  eval(p)
  dev.off()
  print(pltprefix)
}
#' @export
htSortMatrix <- function(M, method="euclidean", ratio=0.5, cutoff=0.25, sort=c(1,2)) {
    if (is.logical(sort)) {
        if (sort) {
            sort = c(1,2)
        } else {
            sort = c()
        }
    }
    rows = 1 %in% sort
    columns = 2 %in% sort
    if (rows) {
        M = order.tsp(M, rows=TRUE, method=method)
    }
    if (columns) {
        M = order.tsp(M, rows=FALSE, method=method)
    }
    if (rows) {
        M = diag.mat3(M, rows=TRUE, ratio=ratio, cutoff=cutoff)
    }
    if (columns) {
        M = diag.mat3(M, rows=FALSE, ratio=ratio, cutoff=cutoff)
    }
    return(M)
}

#' Automatic heatmap
#'
#' @param M matrix to pass
#' @param ux millimeters per grid
#' @param sort Axes to sort upon. Be careful as the matrix is re-sorted before plotting with htSortMatrix
#' @export
autoHeatmap <- function(M, ux=1.5, sort=c(1, 2), method="euclidean", dimname_fontsize=3.5, ratio=0.5, cutoff=0.25, ...) {
    M = htSortMatrix(M, method=method, ratio=ratio, cutoff=cutoff, sort=sort)
    return(ComplexHeatmap::Heatmap(
        M, cluster_rows=FALSE, cluster_columns=FALSE,
        width = ncol(M)*grid::unit(ux, "mm"),
        height = nrow(M)*grid::unit(ux, "mm"),
        row_names_gp=grid::gpar(fontsize=dimname_fontsize),
        column_names_gp=grid::gpar(fontsize=dimname_fontsize),
        ...
    ))
}
#' Save ComplexHeatmap, modified from Carles Boix
#'
#' @param ht ComplexHeatmap heatmap
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param extra Function to call after draw, before figure is saved
#' @param dpi DPI of PNG
#' @export
saveHeatmap <- function(ht, pltprefix, w=7, h=7, dpi=600, extra=NULL, ...) {
    if (w > 100) {
        dpi = 250
    }
    if (!is.null(attr(ht, "mib_w"))) {
        w = attr(ht, "mib_w")
    }
    if (!is.null(attr(ht, "mib_h"))) {
        h = attr(ht, "mib_h")
    }
    pdf(paste0(pltprefix, ".pdf"), width=w, height=h)
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    if (is.function(extra)) {
        extra()
    }
    dev.off()
    png(paste0(pltprefix, ".png"), res=dpi, width=w, height=h, units="in")
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    if (is.function(extra)) {
        extra()
    }
    dev.off()
    print(pltprefix)
}

#' Save ggplot, modified from Carles Boix
#'
#' @param gp GGplot object
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param dpi DPI of PNG
#' @export
saveGGplot <- function(gp, pltprefix, w=7, h=7, dpi=600) {
    ggplot2::ggsave(paste0(pltprefix, ".pdf"), gp, units="in", dpi=dpi, width=w, height=h)
    ggplot2::ggsave(paste0(pltprefix, ".png"), gp, units="in", dpi=dpi, width=w, height=h)
    print(pltprefix)
}

#' Draw triangles in heatmap
#' Pass cell_fun=ht_triangle_split(...) to Heatmap()
#' @param mat.ul Upper left matrix
#' @param mat.lr Lower right matrix
#' @param col.ul Upper left color
#' @param col.lr Lower right color
#' @param lwd lwd for grid::gpar
#' @param ... Arguments to grid::gpar
#' @export
ht_triangle_split <- function(mat.ul, mat.lr, col.ul, col.lr, lwd=0, ...) {
    return(function(j, i, x, y, width, height, fill) {
        ### start with UL, UR, LR, LL
        corners = data.frame(x=as.numeric(width) * c(-1/2, 1/2, 1/2, -1/2),
                             y=as.numeric(height) * c(-1/2, -1/2, 1/2, 1/2),
                             row.names=c("LL", "UL", "UR", "LR"))
        ### upper left
        grid::grid.polygon(x=as.numeric(x) + corners[c("UL", "UR", "LL", "UL"),"x"],
                           y=as.numeric(y) + corners[c("UL", "UR", "LL", "UL"),"y"],
                           gp=grid::gpar(fill=col.ul(mat.ul[i, j]), lwd=lwd, ...))
        grid::grid.polygon(x=as.numeric(x) + corners[c("UR", "LR", "LL", "UR"),"x"],
                           y=as.numeric(y) + corners[c("UR", "LR", "LL", "UR"),"y"],
                           gp=grid::gpar(fill=col.lr(mat.lr[i, j]), lwd=lwd, ...))
    })
}

#' @export
convert_pvalue_to_star <- function(pvalue) {
    stars = rep("", length(pvalue))
    stars[pvalue < 0.05] = "*"
    stars[pvalue < 0.01] = "**"
    stars[pvalue < 0.001] = "***"
    return(setNames(stars, names(pvalue)))
}
#' Draw asterisks
#' Pass cell_fun=ht_asterisks(...) to Heatmap()
#' @param p.matrix P value matrix, unsorted
#' @param gp graphic parameters from grid::gpar(). Recommended to set fontsize
#' @param ... Arguments to grid::grid.text
#' @export
ht_asterisks <- function(p.matrix, gp, ...) {
    return(function(j, i, x, y, w, h, fill) {
        if (is.na(p.matrix[i, j]) | is.nan(p.matrix[i, j])) {
        } else if (p.matrix[i, j] < 0.001) {
            grid::grid.text("***", x, y, gp=gp, ...)
        } else if (p.matrix[i, j] < 0.01) {
            grid::grid.text("**", x, y, gp=gp, ...)
        } else if (p.matrix[i, j] < 0.05) {
            grid::grid.text("*", x, y, gp=gp, ...)
        }
    })
}

#' Draw volcano plot
#' @param df dataframe to pass
#' @param label Text label for points
#' @param title Tile for plot
#' @param threshold.FDR FDR threshold for plotting text
#' @param threshold.log2FC Log2FC threshold for plotting text
#' @param force Force parameter for ggrepel:::geom_text_repel()
#' @param max.overlaps Max overlaps for ggrepel::geom_text_repel()
#' @param quantile.log2FC If less than 1, quantile clip log2FC to fit inside smaller plot area
#' @export
volcano <- function(df, label="gene", title="Volcano plot of",
                    threshold.FDR=0.99999, threshold.log2FC=0.5,
                    force=1, max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                    quantile.log2FC=1) {
    require(ggplot2)
    df$log2FC = qclip(df$log2FC, quantile.log2FC)
    df$FDR = vclip(df$FDR, min=1e-300)
    df$score = -log10(df$FDR) * abs(df$log2FC)
    if (all(label %in% colnames(df))) {
        df$label = apply(df[label], 1, function(x) { paste0(x, collapse=" ") })
    } else {
        df$label = ""
    }
    df$label = ifelse(df$FDR <= threshold.FDR, df$label, rep("", nrow(df)))
    df$label = ifelse(abs(df$log2FC) >= threshold.log2FC, df$label, rep("", nrow(df)))
    g = ggplot(df, aes(x=log2FC, y=-log10(FDR), label=label)) + geom_point() + ggrepel::geom_text_repel(force=force, max.overlaps=max.overlaps) + ggpubr::theme_pubr() + ggtitle(title)
    return(g)
}
