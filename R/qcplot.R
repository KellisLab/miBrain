

#' @export
qcplot <- function(se, assay="counts", n.top.genes=100) {
### 1. Plot top 5 genes
### 2. Plot NG
### 3. Plot total counts
    require(ggplot2)
    M = SummarizedExperiment::assays(se)[[assay]]
    cd = SummarizedExperiment::colData(se)
    cd = as.data.frame(cd)
    cd$n_genes = colSums(M > 0)
    cd$total_counts = colSums(M)
    plt = list(
        "Number of genes"=ggplot(cd, aes(x=n_genes, y=library_id, fill=group)) + scale_fill_brewer(palette="Set2") + geom_bar(stat="identity") + ggpubr::theme_pubr() + ggtitle("Number of genes") + ylab("") + xlab("Number of genes") + scale_y_discrete(labels = NULL, breaks = NULL) + scale_x_continuous(expand=c(0,0)) +  theme(plot.title=element_text(size=11), legend.text=element_text(size=4.5), legend.title=element_text(size=8), axis.ticks.length=unit(.15, "cm")) ,
        "Total counts"=ggplot(cd, aes(x=total_counts, y=library_id, fill=group)) + scale_fill_brewer(palette="Set2") + geom_bar(stat="identity") + ggpubr::theme_pubr() + ggtitle("Total counts") + ylab("") + xlab("Total counts")+ scale_y_discrete(labels = NULL, breaks = NULL) + scale_x_continuous(expand=c(0,0), breaks=scales::pretty_breaks(n=3)) +  theme(plot.title=element_text(size=11), legend.text=element_text(size=4.5), legend.title=element_text(size=8), axis.ticks.length=unit(.25, "cm"))
    )
    
    for (tg in c(100, 500)) {
        tg.idx = apply(M, 2, function(x) {
            head(order(x, decreasing=TRUE), tg)
        })
        cd$top = 100 * sapply(1:nrow(cd), function(i) {
            sum(M[tg.idx[,i],i]) / sum(M[,i])
        })
        cd[[paste0("top", tg)]] = tg
        title = paste0("Percent of counts\nin top ", tg, " genes")
        plt[[title]] = ggplot(cd, aes(x=top, y=library_id, fill=group)) + scale_fill_brewer(palette="Set2") + geom_bar(stat="identity") + ggpubr::theme_pubr() + ggtitle(title) + ylab("") + xlab(paste0("% of counts in top ", tg, " genes")) + scale_y_discrete(labels = NULL, breaks = NULL) + scale_x_continuous(expand=c(0,0)) +  theme(plot.title=element_text(size=11), legend.text=element_text(size=4.5), legend.title=element_text(size=8))
    }
    g = Reduce("+", plt)
    return(g)
}
