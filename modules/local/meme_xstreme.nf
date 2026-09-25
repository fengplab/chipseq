/*
 * De novo + known motif discovery on peak sequences with XSTREME (MEME Suite)
 *
 * All XSTREME defaults are kept except --maxw, which is set per peak set to the length of the shortest
 * read in the final filtered/sorted alignment file(s) that produced the peaks (computed in the workflow).
 */
process MEME_XSTREME {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "bioconda::meme=5.5.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.9--pl5321he99cc7f_1':
        'biocontainers/meme:5.5.9--pl5321he99cc7f_1' }"

    input:
    tuple val(meta), path(fasta), val(maxw)
    path  motif_db

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}_xstreme"), emit: results
    tuple val(meta), path("*.xstreme.html")                       , emit: html, optional: true
    tuple val(meta), path("*.xstreme.combined.meme")              , emit: motifs, optional: true
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def db_args  = motif_db ? [motif_db].flatten().collect { "--m ${it}" }.join(' ') : ''
    """
    if [ ! -s $fasta ]; then
        echo "WARNING: no peak sequences for ${meta.id}; skipping XSTREME" >&2
        mkdir -p ${prefix}_xstreme
    else
        xstreme \\
            --oc ${prefix}_xstreme \\
            --p $fasta \\
            --maxw $maxw \\
            $db_args \\
            $args

        [ -f ${prefix}_xstreme/xstreme.html ]  && cp ${prefix}_xstreme/xstreme.html  ${prefix}.xstreme.html
        [ -f ${prefix}_xstreme/combined.meme ] && cp ${prefix}_xstreme/combined.meme ${prefix}.xstreme.combined.meme
    fi
    printf "maxw_used\\t%s\\n" "${maxw}" > ${prefix}_xstreme/pipeline_maxw.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: \$(xstreme --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_xstreme
    touch ${prefix}.xstreme.html ${prefix}.xstreme.combined.meme
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: 5.5.9
    END_VERSIONS
    """
}
