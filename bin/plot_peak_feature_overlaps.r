#!/usr/bin/env Rscript

################################################
## Summarise/plot per-sample peak overlaps with RepeatMasker / CenSat features
## (features covered >= --min_overlap by a peak)
################################################

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

option_list <- list(
    make_option(c("-i", "--summary_files"), type="character", default=NULL , help="Comma-separated list of *.class_summary.tsv files."),
    make_option(c("-f", "--feature_set"  ), type="character", default="features", help="Feature set name, e.g. 'repeatmasker' or 'censat'."),
    make_option(c("-t", "--top_classes"  ), type="integer"  , default=15   , help="Number of most frequent classes to plot individually."),
    make_option(c("-m", "--min_overlap"  ), type="double"   , default=0.8  , help="Feature coverage threshold used (for labels only)."),
    make_option(c("-o", "--outprefix"    ), type="character", default=NULL , help="Output prefix (default: --feature_set).")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$summary_files)) stop("--summary_files is required")
if (is.null(opt$outprefix)) opt$outprefix <- opt$feature_set

files <- unlist(strsplit(opt$summary_files, ","))
tabs  <- lapply(files, function(f) {
    d <- tryCatch(read.delim(f, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE), error=function(e) NULL)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d
})
tab <- do.call(rbind, tabs[!sapply(tabs, is.null)])
if (is.null(tab) || nrow(tab) == 0) {
    tab <- data.frame(sample=character(), feature_set=character(), class=character(), n_features_covered=integer(),
                      n_peaks=integer(), pct_of_peaks=numeric(), total_peaks=integer())
}
tab <- tab[order(tab$sample, -tab$n_peaks), ]
write.table(tab, paste0(opt$outprefix, ".class_summary.combined.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cls_tab  <- tab[tab$class != "ANY", ]
any_tab  <- tab[tab$class == "ANY", ]
pdf(paste0(opt$outprefix, ".peak_feature_overlap.pdf"), width=11, height=7)
if (nrow(cls_tab) > 0) {
    totals   <- tapply(cls_tab$n_peaks, cls_tab$class, sum)
    top      <- names(sort(totals, decreasing=TRUE))[seq_len(min(opt$top_classes, length(totals)))]
    cls_tab$class_plot <- ifelse(cls_tab$class %in% top, cls_tab$class, "Other")
    agg <- aggregate(cbind(n_peaks, pct_of_peaks) ~ sample + class_plot, data=cls_tab, FUN=sum)
    agg$class_plot <- factor(agg$class_plot, levels=c(top, "Other"))

    thr <- paste0(round(100 * opt$min_overlap), "%")
    p1 <- ggplot(agg, aes(x=sample, y=n_peaks, fill=class_plot)) +
        geom_col(position="dodge") +
        labs(title=paste0(opt$feature_set, ": peaks covering >= ", thr, " of a feature, by class"),
             x=NULL, y="Number of peaks", fill="Class") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
    print(p1)

    p2 <- ggplot(agg, aes(x=sample, y=pct_of_peaks, fill=class_plot)) +
        geom_col(position="dodge") +
        labs(title=paste0(opt$feature_set, ": % of peaks covering >= ", thr, " of a feature, by class"),
             x=NULL, y="% of peaks", fill="Class") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
    print(p2)

    p3 <- ggplot(agg, aes(x=class_plot, y=sample, fill=pct_of_peaks)) +
        geom_tile(colour="white") +
        geom_text(aes(label=n_peaks), size=2.5) +
        scale_fill_gradient(low="white", high="firebrick") +
        labs(title=paste0(opt$feature_set, ": peaks per class (label = count)"), x="Class", y=NULL, fill="% peaks") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
    print(p3)
}
if (nrow(any_tab) > 0) {
    p4 <- ggplot(any_tab, aes(x=sample, y=pct_of_peaks)) +
        geom_col(fill="steelblue") +
        geom_text(aes(label=n_peaks), vjust=-0.3, size=3) +
        labs(title=paste0(opt$feature_set, ": % of peaks covering >= ", round(100 * opt$min_overlap), "% of any feature"),
             x=NULL, y="% of peaks") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
    print(p4)
}
invisible(dev.off())

## MultiQC bargraph (samples x classes, number of peaks)
mqc <- paste0(opt$outprefix, ".class_summary_mqc.tsv")
id  <- paste0("peak_", gsub("[^A-Za-z0-9_]", "_", opt$feature_set), "_overlap")
hdr <- c(
    paste0("# id: '", id, "'"),
    paste0("# section_name: 'MERGED LIB: ", opt$feature_set, " peak annotation'"),
    paste0("# description: 'number of peaks covering at least ", round(100 * opt$min_overlap),
           "% of one or more ", opt$feature_set, " features, split by feature class. A peak can be counted in several classes.'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    paste0("#     id: '", id, "_plot'"),
    paste0("#     title: '", opt$feature_set, ": peaks overlapping features by class'"),
    "#     ylab: 'Number of peaks'"
)
writeLines(hdr, mqc)
if (nrow(cls_tab) > 0) {
    wide <- xtabs(n_peaks ~ sample + class_plot, data=aggregate(n_peaks ~ sample + class_plot, data=cls_tab, FUN=sum))
    wide <- as.data.frame.matrix(wide)
    out  <- cbind(Sample=rownames(wide), wide)
    suppressWarnings(write.table(out, mqc, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE))
} else {
    cat("Sample\tnone\n", file=mqc, append=TRUE)
}
