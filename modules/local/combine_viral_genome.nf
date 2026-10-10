/*
 * Append viral genome(s) to the (renamed) host genome FASTA so reads are aligned to host and viruses together
 */
process COMBINE_VIRAL_GENOME {
    tag "$viral"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path host
    path viral

    output:
    path "*.with_*.fa"       , emit: fasta
    path "viral_contigs.tsv" , emit: contigs
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Name the combined genome after both inputs: the aligner index is named after the FASTA, so an index built
    // with viruses can never be confused with a host-only one
    def host_name  = host.name.replaceAll(/\.gz$/, '').replaceAll(/\.(fa|fasta|fna)$/, '')
    def viral_name = viral.name.replaceAll(/\.gz$/, '').replaceAll(/\.(fa|fasta|fna)$/, '')
    """
    combine_viral_genome.py \\
        --host $host \\
        --viral $viral \\
        --out ${host_name}.with_${viral_name}.fa \\
        --contigs viral_contigs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def host_name  = host.name.replaceAll(/\.gz$/, '').replaceAll(/\.(fa|fasta|fna)$/, '')
    def viral_name = viral.name.replaceAll(/\.gz$/, '').replaceAll(/\.(fa|fasta|fna)$/, '')
    """
    touch ${host_name}.with_${viral_name}.fa viral_contigs.tsv
    echo '"${task.process}": {python: 3.8.3}' > versions.yml
    """
}
