/*
 * Statistically summarise how multi-mapped reads are distributed across repetitive elements.
 *
 * Takes the multi-mapped BAM produced by SAMTOOLS_VIEW_MULTIMAPPED and a BED file of repeat
 * annotations (e.g. a UCSC RepeatMasker `rmsk` track exported to BED, with the repeat
 * name/class in column 4) supplied via `--repeat_masker_bed`, and reports per-interval and
 * per-repeat-class read counts.
 */
process BEDTOOLS_COVERAGE_REPEATS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(bam)
    path  repeat_bed

    output:
    tuple val(meta), path("*.repeat_coverage.txt")      , emit: coverage
    tuple val(meta), path("*.repeat_class_summary.tsv") , emit: summary
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-counts'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Per-interval multi-mapped read counts over each annotated repeat element
    bedtools coverage \\
        $args \\
        -a $repeat_bed \\
        -b $bam \\
        > ${prefix}.repeat_coverage.txt

    # Aggregate to per-repeat-class/family totals (assumes repeat name/class sits in BED col 4,
    # as in a standard RepeatMasker BED export) - number of elements, total and mean
    # multi-mapped read count per class, sorted by total read count descending.
    awk -v OFS='\\t' '
        BEGIN { print "repeat_class", "n_elements", "total_multimapped_reads", "mean_reads_per_element" }
        {
            n[\$4]++
            reads[\$4] += \$NF
        }
        END {
            for (cls in n) {
                printf "%s\\t%d\\t%d\\t%.3f\\n", cls, n[cls], reads[cls], reads[cls] / n[cls]
            }
        }
    ' ${prefix}.repeat_coverage.txt | sort -k3,3nr > ${prefix}.repeat_class_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
