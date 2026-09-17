/*
 * Split a BAM file into uniquely-mapped and multi-mapped reads.
 *
 * The pipeline's normal filtering step (BAMTOOLS_FILTER, see conf/modules.config) discards
 * reads with MAPQ 0 whenever `--keep_multi_map` is false, since bwa/bowtie2 report MAPQ 0 for
 * reads that align equally well to more than one location. This module applies the same MAPQ
 * threshold but keeps the discarded (multi-mapped) reads in their own BAM file, rather than
 * throwing them away, so they can be used downstream (see BEDTOOLS_COVERAGE_REPEATS) to look
 * at how multi-mapped reads distribute across repetitive elements of the genome.
 */
process SAMTOOLS_VIEW_MULTIMAPPED {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.multimapped.bam"), emit: bam
    tuple val(meta), path("*.unique.bam")     , emit: unique_bam
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-q 1' // MAPQ threshold: reads BELOW this are written to -U (multi-mapped)
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        -b \\
        $args \\
        -U ${prefix}.multimapped.bam \\
        -o ${prefix}.unique.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
