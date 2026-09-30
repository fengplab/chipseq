/*
 * Rename chromosome names of a genome FASTA / GTF / BED file with a chromosome alias table
 * (e.g. RefSeq 'NC_060925.1' -> UCSC 'chr1') before any index or downstream file is generated
 */
process RENAME_CHROMS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(input)
    path  alias
    val   column

    output:
    tuple val(meta), path("renamed/*")          , emit: renamed
    tuple val(meta), path("*.rename_report.tsv"), emit: report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def format = meta.id == 'fasta' ? 'fasta' : 'tab'
    def name   = input.name.endsWith('.gz') ? input.name[0..-4] : input.name
    """
    mkdir -p renamed
    rename_chroms.py \\
        $input \\
        $alias \\
        renamed/$name \\
        --format $format \\
        --column '$column' \\
        --report ${meta.id}.rename_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
