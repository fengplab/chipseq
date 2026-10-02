/*
 * Group raw Bowtie2 alignments by read name (required by Allo). Written as compressed BAM to save space:
 * Allo converts its input to SAM internally (in a temporary folder it deletes), so results are identical.
 */
process SAMTOOLS_COLLATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.collate.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # NOTE: do NOT use 'collate -f' (fast mode) here; it discards secondary alignments which Allo needs
    samtools \\
        collate \\
        $args \\
        -@ $task.cpus \\
        -T ./collate_tmp \\
        --output-fmt BAM \\
        -o ${prefix}.collate.bam \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.collate.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
