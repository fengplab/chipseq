/*
 * Convert a RepeatMasker (.out / UCSC rmsk.txt) or CenSat (BED) annotation to a sorted BED6+2 file
 */
process FEATURES_TO_BED {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(annotation)
    path  sizes

    output:
    tuple val(meta), path("*.prepared.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${annotation}" == "${prefix}.prepared.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    features_to_bed.py \\
        $annotation \\
        $sizes \\
        ${prefix}.prepared.bed \\
        --type ${meta.id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
