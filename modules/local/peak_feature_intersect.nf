/*
 * Intersect peaks with RepeatMasker / CenSat features.
 *
 * Outputs (per peak set and feature set):
 *   *.<feature>.all_overlaps.tsv        every peak/feature overlap with bp overlap, % of peak and % of feature covered
 *   *.<feature>.overlap_filtered.tsv    overlaps passing the overlap filter: by default the peak covers >= 80% of the
 *                                       feature OR the feature covers >= 80% of the peak (see 'overlap' input)
 *   *.<feature>.peak_annotation.tsv     one row per peak listing covered features/classes (empty if none)
 *   *.<feature>.class_summary.tsv       per class: covered features, peaks covering >=1 feature, % of all peaks
 */
process PEAK_FEATURE_INTERSECT {
    tag "${meta.id}:${feature.id}"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(peaks), val(feature), path(features)
    val   overlap   // [ feature_min_fraction, peak_min_fraction, mode ('either' | 'both' | 'feature' | 'peak') ]

    output:
    tuple val(meta), val(feature), path("*.all_overlaps.tsv")    , emit: all
    tuple val(meta), val(feature), path("*.overlap_filtered.tsv"), emit: covered
    tuple val(meta), val(feature), path("*.peak_annotation.tsv") , emit: peak_annotation
    tuple val(meta), val(feature), path("*.class_summary.tsv")   , emit: summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = "${task.ext.prefix ?: meta.id}.${feature.id}"
    def (min_feature, min_peak, mode) = overlap
    if (!(mode in ['either', 'both', 'feature', 'peak'])) error "Invalid overlap mode '${mode}': use 'either', 'both', 'feature' or 'peak'"
    """
    # Peaks -> BED4 (chrom, start, end, peak_id); consensus files have a header, MACS3 files may not have names
    awk 'BEGIN { FS = OFS = "\\t" }
        NF >= 3 && \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ {
            id = (NF >= 4 && \$4 != "") ? \$4 : \$1 ":" \$2 "-" \$3
            print \$1, \$2, \$3, id
        }' $peaks \\
        | LC_ALL=C sort -k1,1 -k2,2n > peaks.bed4

    TOTAL_PEAKS=\$(wc -l < peaks.bed4 | tr -d ' ')

    # -wo reports the number of overlapping bp as the last column
    bedtools intersect \\
        -a peaks.bed4 \\
        -b $features \\
        -wo \\
        -sorted \\
        $args \\
        > overlaps.raw

    HEADER="peak_chr\\tpeak_start\\tpeak_end\\tpeak_id\\tpeak_length\\tfeature_chr\\tfeature_start\\tfeature_end\\tfeature_name\\tfeature_score\\tfeature_strand\\tfeature_class\\tfeature_family\\tfeature_length\\toverlap_bp\\tpeak_pct_covered\\tfeature_pct_covered\\tfeature_covered_pass\\tpeak_covered_pass\\tpasses_filter"

    printf "\$HEADER\\n" > ${prefix}.all_overlaps.tsv
    printf "\$HEADER\\n" > ${prefix}.overlap_filtered.tsv

    # feature_covered_pass: the peak covers >= min_feature of the feature
    # peak_covered_pass   : the feature covers >= min_peak of the peak
    # passes_filter       : combined according to mode (either = OR, both = AND, feature / peak = that test only)
    awk -v fmin=${min_feature} -v pmin=${min_peak} -v mode=${mode} -v all=${prefix}.all_overlaps.tsv -v cov=${prefix}.overlap_filtered.tsv '
        BEGIN { FS = OFS = "\\t" }
        {
            plen = \$3 - \$2; flen = \$7 - \$6; ov = \$13
            ppct = (plen > 0) ? 100 * ov / plen : 0
            fpct = (flen > 0) ? 100 * ov / flen : 0
            fpass = (flen > 0 && ov / flen >= fmin)
            ppass = (plen > 0 && ov / plen >= pmin)
            if      (mode == "either")  keep = (fpass || ppass)
            else if (mode == "both")    keep = (fpass && ppass)
            else if (mode == "feature") keep = fpass
            else                        keep = ppass
            pass = keep ? "yes" : "no"
            line = sprintf("%s\\t%s\\t%s\\t%s\\t%d\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%d\\t%d\\t%.2f\\t%.2f\\t%s\\t%s\\t%s", \\
                \$1, \$2, \$3, \$4, plen, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, flen, ov, ppct, fpct, fpass ? "yes" : "no", ppass ? "yes" : "no", pass)
            print line >> all
            if (pass == "yes") print line >> cov
        }' overlaps.raw

    # One row per peak (all peaks kept, features that pass the overlap filter collapsed)
    awk -v feat=${feature.id} '
        BEGIN { FS = OFS = "\\t" }
        FNR == NR {
            if (FNR > 1) {
                k = \$4
                n[k]++
                names[k]   = (k in names)   ? names[k] "," \$9  : \$9
                classes[k] = (k in classes) ? classes[k] "," \$12 : \$12
                if (\$17 + 0 > best[k] + 0) best[k] = \$17
            }
            next
        }
        FNR == 1 { print "peak_id", "peak_chr", "peak_start", "peak_end", feat "_n_covered", feat "_names", feat "_classes", feat "_max_feature_pct_covered" }
        { k = \$4; print k, \$1, \$2, \$3, (k in n) ? n[k] : 0, (k in names) ? names[k] : "NA", (k in classes) ? classes[k] : "NA", (k in best) ? best[k] : "NA" }
        ' ${prefix}.overlap_filtered.tsv peaks.bed4 > ${prefix}.peak_annotation.tsv

    # Class summary
    awk -v total=\$TOTAL_PEAKS -v sample=${meta.id} -v feat=${feature.id} '
        BEGIN { FS = OFS = "\\t" }
        NR > 1 {
            c = \$12
            nf[c]++
            if (!((c SUBSEP \$4) in seen)) { seen[c SUBSEP \$4] = 1; np[c]++ }
            if (!(\$4 in anyp)) { anyp[\$4] = 1; nany++ }
        }
        END {
            print "sample", "feature_set", "class", "n_features_covered", "n_peaks", "pct_of_peaks", "total_peaks"
            for (c in nf) printf "%s\\t%s\\t%s\\t%d\\t%d\\t%.2f\\t%d\\n", sample, feat, c, nf[c], np[c], (total > 0) ? 100 * np[c] / total : 0, total
            printf "%s\\t%s\\t%s\\t%d\\t%d\\t%.2f\\t%d\\n", sample, feat, "ANY", NR - 1, nany, (total > 0) ? 100 * nany / total : 0, total
        }' ${prefix}.overlap_filtered.tsv > ${prefix}.class_summary.tsv

    rm -f overlaps.raw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = "${task.ext.prefix ?: meta.id}.${feature.id}"
    """
    touch ${prefix}.all_overlaps.tsv ${prefix}.overlap_filtered.tsv ${prefix}.peak_annotation.tsv ${prefix}.class_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.30.0
    END_VERSIONS
    """
}
