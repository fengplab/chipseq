/*
 * Extract genomic sequences underlying peak intervals for motif discovery
 */
process PEAKS_TO_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(peaks), val(maxw)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.peaks.fa"), val(maxw), emit: fasta
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Keep chrom/start/end only (consensus files carry a header line and extra columns)
    awk 'BEGIN { FS = OFS = "\\t" } NF >= 3 && \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ && \$3 > \$2 { print \$1, \$2, \$3 }' $peaks \\
        | sort -k1,1 -k2,2n -u > peaks.bed3

    bedtools getfasta \\
        $args \\
        -fi $fasta \\
        -bed peaks.bed3 \\
        -fo ${prefix}.peaks.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.peaks.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.30.0
    END_VERSIONS
    """
}
