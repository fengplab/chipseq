#!/usr/bin/env Rscript

################################################
## Differential binding analysis of consensus peaks with DESeq2 and visualisation of the
## results together with their annotations:
##   - HOMER annotatePeaks (genomic feature category, distance to TSS, nearest gene)
##   - RepeatMasker / CenSat features covered (>= threshold) by each peak
## One run per antibody (consensus peak set). All pairwise contrasts between sample groups
## (sample id without the trailing _REP<n>) are tested.
################################################

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(RColorBrewer)
    library(pheatmap)
})

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , help="featureCounts output for the consensus peaks."),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=7       , help="First column containing sample counts."),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , help="Column containing the interval id."),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , help="Suffix to strip from count column names."),
    make_option(c("-a", "--homer_file"    ), type="character", default=''      , help="Consensus *.boolean.annotatePeaks.txt (optional)."),
    make_option(c("-e", "--feature_files" ), type="character", default=''      , help="Comma-separated *.feature_covered.tsv files (optional)."),
    make_option(c("-n", "--feature_names" ), type="character", default=''      , help="Comma-separated names matching --feature_files."),
    make_option(c("-q", "--fdr"           ), type="double"   , default=0.05    , help="Adjusted p-value threshold."),
    make_option(c("-l", "--lfc"           ), type="double"   , default=0       , help="Absolute log2 fold-change threshold for calling a peak differential."),
    make_option(c("-t", "--top_n"         ), type="integer"  , default=50      , help="Number of top differential peaks in heatmaps."),
    make_option(c("-v", "--vst"           ), type="logical"  , default=TRUE    , help="Use vst (TRUE) or rlog (FALSE) for heatmaps."),
    make_option(c("-p", "--outprefix"     ), type="character", default='differential', help="Output prefix."),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , help="Number of cores.")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$count_file)) stop("--count_file is required")

if (opt$cores > 1 && requireNamespace("BiocParallel", quietly=TRUE)) {
    BiocParallel::register(BiocParallel::MulticoreParam(opt$cores))
    parallel_deseq <- TRUE
} else {
    parallel_deseq <- FALSE
}

dir_colours <- c("Up"="#D7301F", "Down"="#2171B5", "NS"="grey70")

################################################
## Counts and sample groups
################################################

counts_raw <- read.delim(opt$count_file, comment.char="#", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
ids        <- counts_raw[, opt$id_col]
coords     <- data.frame(interval_id=ids, chr=counts_raw[, 2], start=counts_raw[, 3], end=counts_raw[, 4], stringsAsFactors=FALSE)
coords$start <- as.integer(sub(";.*", "", coords$start))   # SAF written from the consensus BED (0-based starts)
coords$end   <- as.integer(sub(";.*", "", coords$end))
coords$chr   <- sub(";.*", "", coords$chr)

counts <- as.matrix(counts_raw[, opt$count_col:ncol(counts_raw), drop=FALSE])
colnames(counts) <- basename(colnames(counts))
if (opt$sample_suffix != '') colnames(counts) <- sub(opt$sample_suffix, "", colnames(counts), fixed=TRUE)
rownames(counts) <- ids
storage.mode(counts) <- "integer"

groups  <- sub("_REP\\d+$", "", colnames(counts))
coldata <- data.frame(row.names=colnames(counts), condition=factor(groups, levels=unique(sort(groups))))
if (nlevels(coldata$condition) < 2) {
    message("Only one sample group found; nothing to compare.")
    quit(save="no", status=0)
}
if (all(table(coldata$condition) < 2)) {
    message("No sample group has replicates; DESeq2 differential testing is not possible.")
    quit(save="no", status=0)
}

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds, parallel=parallel_deseq, quiet=TRUE)
saveRDS(dds, file=paste0(opt$outprefix, ".dds.rds"))
norm_counts <- counts(dds, normalized=TRUE)
write.table(cbind(interval_id=rownames(norm_counts), round(norm_counts, 3)), paste0(opt$outprefix, ".normalised_counts.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
trans <- if (isTRUE(opt$vst) && nrow(dds) >= 1000) vst(dds, blind=FALSE) else if (isTRUE(opt$vst)) varianceStabilizingTransformation(dds, blind=FALSE) else rlog(dds, blind=FALSE)
trans_mat <- assay(trans)

################################################
## Annotations
################################################

annot <- coords

## HOMER
has_homer <- FALSE
if (opt$homer_file != '' && file.exists(opt$homer_file) && file.info(opt$homer_file)$size > 0) {
    homer <- read.delim(opt$homer_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote="")
    idc   <- if ("interval_id" %in% colnames(homer)) "interval_id" else colnames(homer)[1]
    getcol <- function(nm) if (nm %in% colnames(homer)) homer[[nm]] else rep(NA, nrow(homer))
    h <- data.frame(interval_id=homer[[idc]],
                    homer_annotation=getcol("Annotation"),
                    homer_distance_to_tss=suppressWarnings(as.numeric(getcol("Distance to TSS"))),
                    homer_gene_name=getcol("Gene Name"),
                    stringsAsFactors=FALSE)
    h$homer_category <- trimws(sub("\\s*\\(.*$", "", h$homer_annotation))
    h$homer_category[is.na(h$homer_category) | h$homer_category == ""] <- "NA"
    annot <- merge(annot, h, by="interval_id", all.x=TRUE, sort=FALSE)
    annot$homer_category[is.na(annot$homer_category)] <- "NA"
    has_homer <- TRUE
}

## RepeatMasker / CenSat (features covered >= threshold by a peak)
feature_files <- if (opt$feature_files != '') unlist(strsplit(opt$feature_files, ",")) else character(0)
feature_names <- if (opt$feature_names != '') unlist(strsplit(opt$feature_names, ",")) else character(0)
if (length(feature_names) != length(feature_files)) feature_names <- sub("\\..*$", "", basename(feature_files))
feature_long <- list()
for (i in seq_along(feature_files)) {
    f  <- feature_files[i]
    nm <- feature_names[i]
    d  <- tryCatch(read.delim(f, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote=""), error=function(e) NULL)
    if (is.null(d)) next
    if (nrow(d) > 0) {
        fl <- unique(data.frame(interval_id=d$peak_id, class=d$feature_class, stringsAsFactors=FALSE))
        feature_long[[nm]] <- fl
        coll <- aggregate(class ~ interval_id, data=fl, FUN=function(x) paste(sort(unique(x)), collapse=","))
        nms  <- aggregate(feature_name ~ peak_id, data=d, FUN=function(x) paste(unique(x), collapse=","))
        colnames(nms) <- c("interval_id", paste0(nm, "_features"))
        colnames(coll) <- c("interval_id", paste0(nm, "_classes"))
        annot <- merge(annot, coll, by="interval_id", all.x=TRUE, sort=FALSE)
        annot <- merge(annot, nms,  by="interval_id", all.x=TRUE, sort=FALSE)
    } else {
        feature_long[[nm]] <- data.frame(interval_id=character(), class=character())
        annot[[paste0(nm, "_classes")]]  <- NA
        annot[[paste0(nm, "_features")]] <- NA
    }
}
rownames(annot) <- annot$interval_id
annot <- annot[rownames(dds), ]

################################################
## Helper functions
################################################

## Fisher test of class membership in Up/Down peaks against all tested peaks
enrichment <- function(member_ids_by_class, res_df) {
    out <- list()
    universe <- res_df$interval_id
    for (dirn in c("Up", "Down")) {
        in_dir <- res_df$interval_id[res_df$direction == dirn]
        if (length(in_dir) == 0) next
        for (cl in names(member_ids_by_class)) {
            mem <- intersect(member_ids_by_class[[cl]], universe)
            a <- sum(in_dir %in% mem); b <- length(in_dir) - a
            c <- length(mem) - a;      d <- length(universe) - a - b - c
            ft <- fisher.test(matrix(c(a, b, c, d), nrow=2))
            out[[length(out) + 1]] <- data.frame(direction=dirn, class=cl, n_in_direction=a, n_direction=length(in_dir),
                                                 n_background=length(mem), n_universe=length(universe),
                                                 pct_in_direction=100 * a / length(in_dir), pct_background=100 * length(mem) / length(universe),
                                                 odds_ratio=unname(ft$estimate), pvalue=ft$p.value, stringsAsFactors=FALSE)
        }
    }
    if (length(out) == 0) return(NULL)
    e <- do.call(rbind, out)
    e$padj <- p.adjust(e$pvalue, method="BH")
    e
}

plot_enrichment <- function(e, title) {
    if (is.null(e) || nrow(e) == 0) return(invisible(NULL))
    e$log2OR <- log2(pmin(pmax(e$odds_ratio, 1 / 64), 64))
    e$sig    <- ifelse(e$padj < 0.05, "padj < 0.05", "n.s.")
    ggplot(e, aes(x=reorder(class, log2OR), y=log2OR, fill=direction, alpha=sig)) +
        geom_col(position=position_dodge(width=0.8), width=0.75) +
        geom_hline(yintercept=0, size=0.3) +
        scale_fill_manual(values=dir_colours) +
        scale_alpha_manual(values=c("padj < 0.05"=1, "n.s."=0.35)) +
        coord_flip() +
        labs(title=title, x=NULL, y="log2 odds ratio vs all consensus peaks (Fisher's exact test)", alpha=NULL, fill=NULL) +
        theme_bw()
}

composition_plot <- function(df, cat_col, title, top=15) {
    tab <- as.data.frame(table(direction=df$direction, category=df[[cat_col]]), stringsAsFactors=FALSE)
    tab <- tab[tab$Freq > 0, ]
    if (nrow(tab) == 0) return(invisible(NULL))
    tops <- names(sort(tapply(tab$Freq, tab$category, sum), decreasing=TRUE))[seq_len(min(top, length(unique(tab$category))))]
    tab$category <- ifelse(tab$category %in% tops, tab$category, "Other")
    tab <- aggregate(Freq ~ direction + category, data=tab, FUN=sum)
    totals <- table(df$direction)
    tab$direction <- factor(paste0(tab$direction, "\n(n=", totals[tab$direction], ")"),
                            levels=paste0(names(dir_colours), "\n(n=", totals[names(dir_colours)], ")"))
    ncat <- length(unique(tab$category))
    pal  <- if (ncat <= 12) brewer.pal(max(3, ncat), "Set3")[seq_len(ncat)] else colorRampPalette(brewer.pal(12, "Set3"))(ncat)
    ggplot(tab, aes(x=direction, y=Freq, fill=category)) +
        geom_col(position="fill") +
        scale_y_continuous(labels=function(x) paste0(100 * x, "%")) +
        scale_fill_manual(values=pal) +
        labs(title=title, x=NULL, y="Proportion of peaks", fill=NULL) +
        theme_bw()
}

################################################
## Pairwise contrasts
################################################

lvls      <- levels(coldata$condition)
pairs     <- combn(lvls, 2, simplify=FALSE)
summary_l <- list()
ma_by_contrast <- list()

for (pr in pairs) {
    g1 <- pr[1]; g2 <- pr[2]
    contrast_name <- paste0(g2, "vs", g1)
    cprefix <- paste0(opt$outprefix, ".", contrast_name)

    res <- results(dds, contrast=c("condition", g2, g1), alpha=opt$fdr, lfcThreshold=0, parallel=parallel_deseq)
    res_df <- data.frame(interval_id=rownames(res), as.data.frame(res), stringsAsFactors=FALSE)
    res_df$direction <- "NS"
    sig <- !is.na(res_df$padj) & res_df$padj < opt$fdr & abs(res_df$log2FoldChange) >= opt$lfc
    res_df$direction[sig & res_df$log2FoldChange > 0] <- "Up"
    res_df$direction[sig & res_df$log2FoldChange < 0] <- "Down"
    res_df$direction <- factor(res_df$direction, levels=names(dir_colours))
    res_df <- cbind(res_df, annot[res_df$interval_id, setdiff(colnames(annot), "interval_id"), drop=FALSE])
    res_df <- cbind(res_df, round(norm_counts[res_df$interval_id, , drop=FALSE], 3))
    res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last=TRUE), ]

    write.table(res_df, paste0(cprefix, ".results.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(res_df[res_df$direction != "NS", ], paste0(cprefix, ".significant.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    sig_bed <- res_df[res_df$direction != "NS", c("chr", "start", "end", "interval_id", "log2FoldChange", "direction")]
    write.table(sig_bed[order(sig_bed$chr, sig_bed$start), ], paste0(cprefix, ".significant.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    n_up   <- sum(res_df$direction == "Up");   n_down <- sum(res_df$direction == "Down")
    summary_l[[contrast_name]] <- data.frame(contrast=contrast_name, group_numerator=g2, group_denominator=g1,
                                             n_tested=nrow(res_df), n_up=n_up, n_down=n_down, fdr=opt$fdr, lfc=opt$lfc)
    subtitle <- paste0(g2, " vs ", g1, " | padj < ", opt$fdr, if (opt$lfc > 0) paste0(", |log2FC| >= ", opt$lfc) else "",
                       " | Up: ", n_up, "  Down: ", n_down)

    pdf(paste0(cprefix, ".plots.pdf"), width=10, height=8)

    ## 1. MA plot
    ma <- res_df[!is.na(res_df$log2FoldChange), ]
    print(ggplot(ma, aes(x=log10(baseMean + 1), y=log2FoldChange, colour=direction)) +
        geom_point(size=0.7, alpha=0.6) + geom_hline(yintercept=0, size=0.3) +
        scale_colour_manual(values=dir_colours) +
        labs(title=paste0("MA plot: ", contrast_name), subtitle=subtitle, x="log10(mean normalised count + 1)", y="log2 fold change", colour=NULL) +
        theme_bw())

    ## 2. Volcano plot
    vol <- res_df[!is.na(res_df$padj), ]
    vol$mlog10 <- -log10(pmax(vol$padj, .Machine$double.xmin))
    top_lab <- head(vol[vol$direction != "NS", ], 15)
    lab_col <- if (has_homer) "homer_gene_name" else "interval_id"
    p_vol <- ggplot(vol, aes(x=log2FoldChange, y=mlog10, colour=direction)) +
        geom_point(size=0.8, alpha=0.6) +
        geom_hline(yintercept=-log10(opt$fdr), linetype="dashed", size=0.3) +
        scale_colour_manual(values=dir_colours) +
        labs(title=paste0("Volcano plot: ", contrast_name), subtitle=subtitle, x="log2 fold change", y="-log10 adjusted p-value", colour=NULL) +
        theme_bw()
    if (opt$lfc > 0) p_vol <- p_vol + geom_vline(xintercept=c(-opt$lfc, opt$lfc), linetype="dashed", size=0.3)
    if (nrow(top_lab) > 0) p_vol <- p_vol + geom_text(data=top_lab, aes(label=.data[[lab_col]]), size=2.5, vjust=-0.6, show.legend=FALSE, check_overlap=TRUE)
    print(p_vol)

    ## 3. Heatmap of the top differential peaks
    top_ids <- head(res_df$interval_id[res_df$direction != "NS"], opt$top_n)
    if (length(top_ids) >= 2) {
        samples <- rownames(coldata)[coldata$condition %in% c(g1, g2)]
        mat <- trans_mat[top_ids, samples, drop=FALSE]
        mat <- t(scale(t(mat))); mat[is.na(mat)] <- 0
        row_ann <- data.frame(direction=as.character(res_df[top_ids, "direction"]), row.names=top_ids)
        if (has_homer) row_ann$HOMER <- annot[top_ids, "homer_category"]
        for (nm in names(feature_long)) {
            cc <- annot[top_ids, paste0(nm, "_classes")]
            row_ann[[nm]] <- ifelse(is.na(cc), "none", sub(",.*", "", cc))
        }
        col_ann <- data.frame(condition=as.character(coldata[samples, "condition"]), row.names=samples)
        rn <- if (has_homer) paste0(top_ids, " (", ifelse(is.na(annot[top_ids, "homer_gene_name"]), "", annot[top_ids, "homer_gene_name"]), ")") else top_ids
        pheatmap(mat, annotation_row=row_ann, annotation_col=col_ann, labels_row=rn, fontsize_row=6,
                 color=colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), cluster_cols=TRUE,
                 main=paste0("Top ", length(top_ids), " differential peaks (row z-score): ", contrast_name))
    }

    ## 4. HOMER annotation of differential peaks
    if (has_homer) {
        print(composition_plot(res_df, "homer_category", paste0("HOMER genomic annotation by direction: ", contrast_name)))
        cats <- split(res_df$interval_id, res_df$homer_category)
        e_h  <- enrichment(cats, res_df)
        if (!is.null(e_h)) {
            write.table(e_h, paste0(cprefix, ".homer_enrichment.txt"), sep="\t", quote=FALSE, row.names=FALSE)
            print(plot_enrichment(e_h, paste0("HOMER category enrichment in differential peaks: ", contrast_name)))
        }
        dd <- res_df[!is.na(res_df$homer_distance_to_tss), ]
        if (nrow(dd) > 0) {
            print(ggplot(dd, aes(x=sign(homer_distance_to_tss) * log10(abs(homer_distance_to_tss) + 1), colour=direction)) +
                geom_density(size=0.8) + scale_colour_manual(values=dir_colours) +
                labs(title=paste0("Distance to nearest TSS (HOMER): ", contrast_name), x="signed log10(distance to TSS + 1)", y="Density", colour=NULL) +
                theme_bw())
        }
    }

    ## 5. RepeatMasker / CenSat annotation of differential peaks
    for (nm in names(feature_long)) {
        fl <- feature_long[[nm]]
        fl <- fl[fl$interval_id %in% res_df$interval_id, ]
        overlap_any <- res_df$interval_id %in% fl$interval_id
        any_df <- data.frame(direction=res_df$direction, covered=ifelse(overlap_any, "covers >=1 feature", "none"))
        any_tab <- as.data.frame(prop.table(table(any_df$direction, any_df$covered), 1))
        colnames(any_tab) <- c("direction", "covered", "prop")
        print(ggplot(any_tab, aes(x=direction, y=prop, fill=covered)) + geom_col() +
            scale_fill_manual(values=c("covers >=1 feature"="darkorange", "none"="grey80")) +
            scale_y_continuous(labels=function(x) paste0(100 * x, "%")) +
            labs(title=paste0(nm, ": proportion of peaks covering >= threshold of a feature: ", contrast_name), x=NULL, y="Proportion of peaks", fill=NULL) +
            theme_bw())

        if (nrow(fl) > 0) {
            long <- merge(fl, res_df[, c("interval_id", "direction")], by="interval_id")
            per_dir <- aggregate(interval_id ~ direction + class, data=long, FUN=function(x) length(unique(x)))
            colnames(per_dir)[3] <- "n_peaks"
            tot <- table(res_df$direction)
            per_dir$pct <- 100 * per_dir$n_peaks / as.numeric(tot[as.character(per_dir$direction)])
            tops <- names(sort(tapply(per_dir$n_peaks, per_dir$class, sum), decreasing=TRUE))[seq_len(min(20, length(unique(per_dir$class))))]
            per_dir <- per_dir[per_dir$class %in% tops, ]
            per_dir$class <- factor(per_dir$class, levels=rev(tops))
            print(ggplot(per_dir, aes(x=class, y=pct, fill=direction)) +
                geom_col(position=position_dodge(width=0.85), width=0.8) + coord_flip() +
                scale_fill_manual(values=dir_colours) +
                labs(title=paste0(nm, " classes covered by peaks, by direction: ", contrast_name),
                     x=NULL, y="% of peaks in direction", fill=NULL) + theme_bw())

            members <- split(fl$interval_id, fl$class)
            e_f <- enrichment(members, res_df)
            if (!is.null(e_f)) {
                write.table(e_f, paste0(cprefix, ".", nm, "_enrichment.txt"), sep="\t", quote=FALSE, row.names=FALSE)
                print(plot_enrichment(e_f[e_f$class %in% tops, ], paste0(nm, " class enrichment in differential peaks: ", contrast_name)))
            }

            vol$feature <- ifelse(vol$interval_id %in% fl$interval_id, paste0("covers ", nm, " feature"), "no feature")
            print(ggplot(vol, aes(x=log2FoldChange, y=mlog10)) +
                geom_point(data=vol[vol$feature == "no feature", ], colour="grey80", size=0.6) +
                geom_point(data=vol[vol$feature != "no feature", ], aes(colour=direction), size=0.9, alpha=0.8) +
                geom_hline(yintercept=-log10(opt$fdr), linetype="dashed", size=0.3) +
                scale_colour_manual(values=dir_colours) +
                labs(title=paste0("Volcano plot highlighting peaks covering ", nm, " features: ", contrast_name),
                     subtitle="grey = no covered feature", x="log2 fold change", y="-log10 adjusted p-value", colour=NULL) +
                theme_bw())
        }
    }
    invisible(dev.off())
}

################################################
## Summary across contrasts (+ MultiQC)
################################################

summary_df <- do.call(rbind, summary_l)
write.table(summary_df, paste0(opt$outprefix, ".differential_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE)

pdf(paste0(opt$outprefix, ".differential_summary.pdf"), width=9, height=6)
sl <- rbind(data.frame(contrast=summary_df$contrast, direction="Up",   n=summary_df$n_up),
            data.frame(contrast=summary_df$contrast, direction="Down", n=summary_df$n_down))
print(ggplot(sl, aes(x=contrast, y=n, fill=direction)) + geom_col(position="dodge") +
    geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.3, size=3) +
    scale_fill_manual(values=dir_colours) +
    labs(title=paste0(opt$outprefix, ": differential consensus peaks (padj < ", opt$fdr, ")"), x=NULL, y="Number of peaks", fill=NULL) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)))
invisible(dev.off())

mqc <- paste0(opt$outprefix, ".differential_mqc.tsv")
id  <- paste0("differential_peaks_", gsub("[^A-Za-z0-9_]", "_", opt$outprefix))
writeLines(c(
    paste0("# id: '", id, "'"),
    paste0("# section_name: 'MERGED LIB: ", opt$outprefix, " DESeq2 differential peaks'"),
    paste0("# description: 'number of consensus peaks with significantly higher (Up) or lower (Down) signal in the first group of each contrast (padj < ",
           opt$fdr, if (opt$lfc > 0) paste0(", |log2FC| >= ", opt$lfc) else "", ").'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    paste0("#     id: '", id, "_plot'"),
    paste0("#     title: '", opt$outprefix, ": differential consensus peaks'"),
    "#     ylab: 'Number of peaks'",
    "# categories:",
    "#     Up:",
    "#         color: '#D7301F'",
    "#     Down:",
    "#         color: '#2171B5'",
    "Contrast\tUp\tDown"), mqc)
write.table(summary_df[, c("contrast", "n_up", "n_down")], mqc, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

message("Done.")
