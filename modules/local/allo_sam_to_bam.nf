/*
 * Post-process Allo output:
 *   - Allo keeps the SAM record of the chosen alignment untouched, so an allocated multi-mapper that
 *     Bowtie2 reported as a *secondary* alignment (FLAG 0x100) would still be flagged secondary and would
 *     be silently dropped by Picard/MACS3/featureCounts. The 0x100 bit is cleared for records carrying the
 *     Allo ZA/ZZ tags so that exactly one primary record per read (pair) is carried forward.
 *   - Convert to BAM; coordinate sorting/indexing/stats are done by BAM_SORT_STATS_SAMTOOLS afterwards.
 */
process ALLO_SAM_TO_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(allo_bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${allo_bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools view -h $allo_bam \\
        | awk 'BEGIN { FS = OFS = "\\t" }
        /^@/ { print; next }
        {
            if (\$0 ~ /\\tZ[AZ]:Z:[0-9]+/ && int(\$2 / 256) % 2 == 1) { \$2 = \$2 - 256 }
            print
        }' \\
        | samtools view $args -@ $task.cpus -b -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
