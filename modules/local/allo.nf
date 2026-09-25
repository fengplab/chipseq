/*
 * Allocate multi-mapped reads with Allo (https://github.com/seqcode/allo)
 *
 * Input : raw, read-name-grouped SAM produced by Bowtie2 (-k N) + samtools collate
 * Output: Allo SAM (one alignment per read / pair, allocated multi-mappers tagged ZA/ZZ)
 *
 * NOTE: the 'container' and 'conda' directives for this process are intentionally NOT declared here.
 *       They are set in conf/modules.config so that '--allo_use_local' can remove them entirely and
 *       force Nextflow to run the task with the 'allo' executable found on the host $PATH.
 *       Default container: biocontainers/allo:1.2.0--pyhdfd78af_0 (bioconda::allo=1.2.0)
 */
process ALLO {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.allo.sam")   , emit: sam
    tuple val(meta), path("*.allo.log")   , emit: log
    tuple val(meta), path("*.allo_mqc.tsv"), emit: mqc
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq    = meta.single_end ? 'se' : 'pe'
    def unit   = meta.single_end ? 'reads' : 'pairs'
    if ("${sam}" == "${prefix}.allo.sam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    if ! command -v allo >/dev/null 2>&1; then
        echo "ERROR: 'allo' executable not found. Either run with a container/conda profile or install Allo locally (pip install bio-allo) when using --allo_use_local." >&2
        exit 1
    fi

    # TensorFlow is chatty and may try to use every core on the node; keep it within the task allocation
    export TF_CPP_MIN_LOG_LEVEL=2
    export OMP_NUM_THREADS=1
    export TF_NUM_INTRAOP_THREADS=1
    export TF_NUM_INTEROP_THREADS=1

    allo \\
        $sam \\
        -seq $seq \\
        -p $task.cpus \\
        -o ${prefix}.allo.sam \\
        $args \\
        2>&1 | tee ${prefix}.allo.log

    # Allo exits 0 on some fatal errors; make sure it actually produced output
    if [ ! -s ${prefix}.allo.sam ] || ! grep -q "Allocation finished" ${prefix}.allo.log; then
        echo "ERROR: Allo did not finish successfully, see ${prefix}.allo.log" >&2
        exit 1
    fi

    # MultiQC custom content (one row per library, headers are added when the files are collated)
    awk -v s="${meta.id}" -v unit="$unit" '
        /Total uniquely mapped/              { split(\$0,a,": "); um=a[2] }
        /Total number of (reads|pairs) allocated/ { split(\$0,a,": "); al=a[2] }
        /Total number of reads filtered/     { split(\$0,a,": "); fi=a[2] }
        /Average number of alignments/       { split(\$0,a,": "); av=a[2] }
        END {
            tot = um + al
            pct = (tot > 0) ? sprintf("%.2f", 100 * al / tot) : "0"
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", s, unit, um, al, pct, fi, av
        }' ${prefix}.allo.log > ${prefix}.allo_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allo: \$(allo --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.allo.sam ${prefix}.allo.log
    printf "${meta.id}\\treads\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}.allo_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allo: 1.2.0
    END_VERSIONS
    """
}
