/*
 * Build a reusable BigBed of RepeatMasker / CenSat features (pass it back in with --repeatmasker_bigbed / --censat_bigbed)
 */
process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ucsc-bedtobigbed=482"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0' :
        'biocontainers/ucsc-bedtobigbed:482--hdc0a859_0' }"

    input:
    tuple val(meta), path(bed)
    path  sizes
    path  autosql

    output:
    tuple val(meta), path("*.bb"), emit: bigbed
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '482' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bedToBigBed \\
        -type=bed6+2 \\
        -as=$autosql \\
        -tab \\
        $args \\
        $bed \\
        $sizes \\
        ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
