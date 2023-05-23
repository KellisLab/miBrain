
#' Load STAR ReadsPerGene into SummarizedExperiment
#'
#' @param star.prefix Vector of directories with STAR ReadsPerGene.out.tab inside. Can be named
#' @param index STAR index used. Will use geneInfo.tab to load in gene names
#' @param strand By default, use unstranded. Can use either first or second strand as well.
#' @param featureCounts use featureCounts output file if exists. Will use instead of STAR ReadsPerGene.out.tab
#' @export
read_star <- function(star.prefix, index="/net/bmc-lab5/data/kellis/group/Benjamin/ref/STAR_gencode43/", featureCounts="featureCounts", strand=NULL, ...) {
    gf = NULL
    if (file.exists(paste0(index, "/geneInfo.tab"))) {
        gf = read.table(paste0(index, "/geneInfo.tab"), skip=1,sep="\t")
        if (ncol(gf) == 3) {
            colnames(gf) = c("gene_id", "gene_name", "gene_type")
        } else {
            colnames(gf) = "gene_id"
        }
        rownames(gf) = gf$gene_id
    }
    xf = NULL
    for (i in seq_along(star.prefix)) {
        if (all(file.exists(paste0(star.prefix, featureCounts)))) {
            print(paste0("Using ", star.prefix[i], "featureCounts"))
            ### Only use featureCounts if exists for all in prefix.
            df = read.table(paste0(star.prefix[i], featureCounts), sep="\t", skip="#", header=TRUE)
            if (ncol(df) == 7) {
                df = df[c("Geneid", colnames(df)[ncol(df)])]
                colnames(df) = c("gene_id", "unstranded")
            } else {
                stop(paste0("Incorrect number of columns in ", star.prefix[i], featureCounts, ", ", ncol(df)))
            }
        } else {
            print(paste0("Using ", star.prefix[i], "ReadsPerGene.out.tab"))
            df = read.table(paste0(star.prefix[i], "ReadsPerGene.out.tab"), sep="\t", col.names=c("gene_id", "unstranded", "first", "second"))
        }
        df$sample = ifelse(is.null(names(star.prefix)), basename(star.prefix[i]), names(star.prefix)[i])
        if (is.null(xf)) {
            xf = df
        } else {
            xf = rbind(xf, df)
        }
    }
    if (is.null(strand)) {
        xf$count = xf$unstranded
    } else {
        xf$count = xf[[strand]]
    }
    if (is.null(gf)) {
        M = pivot(xf, "gene_id", "sample", "count")
        se = SummarizedExperiment::SummarizedExperiment(list(counts=M), ...)
    } else {
        M = pivot(xf[xf$gene_id %in% gf$gene_id,], "gene_id", "sample", "count")
        se = SummarizedExperiment::SummarizedExperiment(list(counts=M), rowData=gf[rownames(M),], ...)
        if ("gene_name" %in% colnames(gf)) {
            rownames(se) = var_names_make_unique(gf[rownames(se),]$gene_name)
        }
    }
    return(se)
}

#' List BMC filenames in batch directory
#' @param df Dataframe with rownames for sample names, and colnames for colData columns
#' @param seq.dir Outer directory where sample names are located
#' @param extra Boolean indicating whether BMC style suffixes should be added
#' @export
bulk_aggregate_star <- function(df, seq.dir, extra=TRUE, ...) {
    sample.dir.list = sapply(rownames(df), function(rn) {
        sample.dir = list.files(paste0(seq.dir, "/"),
                                pattern=paste0("^", rn,
                                               ifelse(extra, "-[A-Z0-9]+", ""),
                                               "$"),
                                full.names=TRUE)
        stopifnot(length(sample.dir) == 1)
        return(sample.dir)
    })
    stopifnot(length(sample.dir.list) == nrow(df))
    se = read_star(setNames(paste0(sample.dir.list, "/"), rownames(df)), ...)
    for (i in seq_along(colnames(df))) {
        SummarizedExperiment::colData(se)[[ colnames(df)[i] ]] = df[[ i ]]
    }
    return(se)
}


#' Add genomicranges of gene to SummarizedExperiment
#' @param se SummarizedExperiment object
#' @param gtf GTF/GFF file containing positions
#' @param by Gene id to split on
#' @export
se_gene_ranges <- function(se, gtf, by="gene_name") {
    rd = SummarizedExperiment::rowData(se)
    gf = as.data.frame(rtracklayer::readGFF(gtf))
    gf = gf[gf$type == "gene",]
    rownames(gf) = var_names_make_unique(gf[[by]])
    gr = with(gf, GenomicRanges::GRanges(seqid, ranges=IRanges::IRanges(start, end)))
    names(gr) = rownames(gf)
    gr = gr[var_names_make_unique(rd[[by]])]
    for (name in colnames(rd)) {
        S4Vectors::mcols(gr)[[name]] = rd[[name]]
    }
    names(gr) = rownames(rd)
    SummarizedExperiment::rowRanges(se) = gr
    rd$gene_length = GenomicRanges::width(gr)
    SummarizedExperiment::rowData(se) = rd
    return(se)
}

#' Count (from featureCounts) raw data from a spreadsheet into SummarizedExperiment object
#' @param samplesheet.csv CSV with rownames in 1st column to be used as colData
#' @param seq.dir Directory where batch subdir lives
#' @param star.index STAR index passed to read_star
#' @param gtf GTF/GFF file to add gene_length column
#' @return SummarizedExperiment object across all batches
#' @export
bulk_aggregate_counts <- function(samplesheet.csv, seq.dir, star.index="/net/bmc-lab5/data/kellis/group/Benjamin/ref/STAR_gencode43/", gtf="/net/bmc-lab5/data/kellis/group/Benjamin/ref/gencode.v43.annotation.gff3.gz") {
    df = read.csv(samplesheet.csv, row.names=1)
    df$library_id = rownames(df)
    all.se = do.call(SummarizedExperiment::cbind, lapply(split(df, df$batch), function(bf) {
        batch.dir = paste0(seq.dir, "/", bf$batch[1], "/")
        return(bulk_aggregate_star(bf, batch.dir, index=star.index))
    }))
    colnames(all.se) = make.unique(all.se$title)
    if (file.exists(as.character(gtf))) {
        all.se = se_gene_ranges(all.se, gtf=gtf)
        SummarizedExperiment::assays(all.se)[["RPKM"]] = edgeR::rpkm(all.se, gene.length="gene_length")
    }
    return(all.se)
}

#' @export
bulk_aggregate_comparison <- function(se, comparison.csv, logTMM=FALSE) {
    df = read.csv(comparison.csv)
    return(lapply(split(df, df$comparison), function(cf) {
        case = cf[cf$case == TRUE, c("group", "batch")]
        stopifnot(length(unique(case$group))==1)
        ctrl = cf[cf$case == FALSE, c("group", "batch")]
        stopifnot(length(unique(ctrl$group))==1)
        se.case = se[,(se$batch %in% case$batch) & (se$group %in% case$group)]
        se.ctrl = se[,(se$batch %in% ctrl$batch) & (se$group %in% ctrl$group)]
        cc = SummarizedExperiment::cbind(se.case, se.ctrl)
        SummarizedExperiment::metadata(cc) = list(comparison=unique(cf$comparison),
                                                  case=unique(case$group),
                                                  control=unique(ctrl$group))
        return(se_tmm(cc, log=logTMM))
    }))
}
