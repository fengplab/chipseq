/*
 * Allocate multi-mapped reads with Allo (https://github.com/seqcode/allo)
 *
 * Input : raw, read-name-grouped Bowtie2 (-k N) alignments from samtools collate (BAM; Allo converts it to SAM
 *         internally in a temporary folder that it deletes)
 * Output: Allo alignments (one per read / pair, allocated multi-mappers tagged ZA/ZZ), compressed to BAM
 *         inside the task so the large uncompressed SAM does not stay in the work directory
 *
 * NOTE: the 'container' and 'conda' directives for this process are intentionally NOT declared here.
 *       They are set in conf/modules.config so that '--allo_use_local' can remove them entirely and
 *       force Nextflow to run the task with the 'allo' executable found on the host $PATH.
 *       Container: build from containers/allo/ (the public biocontainer lacks 'keras' and cannot import Allo).
 */
process ALLO {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(collated)

    output:
    tuple val(meta), path("*.allo.bam")   , emit: bam
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
    // Every Allo worker loads TensorFlow + the CNN and receives a full copy of the genome-wide read landscape,
    // so memory grows with the worker count. Cap it (--allo_max_workers) and halve it on each retry.
    def max_workers = params.allo_max_workers ? (params.allo_max_workers as int) : task.cpus
    def workers     = Math.max(1, (Math.min(task.cpus as int, max_workers) / (1 << (task.attempt - 1))) as int)
    if ("${collated}" == "${prefix}.allo.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    if ! command -v allo >/dev/null 2>&1; then
        echo "ERROR: 'allo' executable not found. Either pass --allo_container / use -profile conda, or install Allo locally (pip install bio-allo keras) when using --allo_use_local." >&2
        exit 1
    fi

    # TensorFlow is chatty and may try to use every core on the node; keep it within the task allocation
    export TF_CPP_MIN_LOG_LEVEL=2
    export OMP_NUM_THREADS=1
    export TF_NUM_INTRAOP_THREADS=1
    export TF_NUM_INTEROP_THREADS=1

    # Work around a segfault in TensorFlow's SavedModel fingerprinting (CreateFingerprintDef), hit while Allo
    # converts its CNN to TF Lite in every worker, seen with pip TensorFlow inside conda envs (protobuf mismatch).
    # The fingerprint is only metadata of a temporary SavedModel, so disabling it does not change Allo's results.
    # sitecustomize.py is executed by every Python process of this task, including Allo's joblib workers.
    if [ "${params.allo_tf_fingerprint_workaround}" == "true" ]; then
        mkdir -p .allo_tf_workaround
        cat > .allo_tf_workaround/sitecustomize.py <<-'END_PY'
    try:
        from tensorflow.core.config import flags as _tf_flags
        _tf_flags.config().saved_model_fingerprinting.reset(False)
    except Exception:
        pass
    END_PY
        export PYTHONPATH="\$PWD/.allo_tf_workaround\${PYTHONPATH:+:\$PYTHONPATH}"
    fi

    set +e
    allo \\
        $collated \\
        -seq $seq \\
        -p $workers \\
        -o ${prefix}.allo.sam \\
        $args \\
        2>&1 | tee ${prefix}.allo.log
    ALLO_RC=\${PIPESTATUS[0]}
    set -e

    if grep -qE "No module named '(tensorflow[.]keras|keras)'" ${prefix}.allo.log; then
        echo "ERROR: the Allo environment is missing the 'keras' package required by TensorFlow >= 2.16." >&2
        echo "       Install it ('pip install keras') or build the container from containers/allo/ and pass --allo_container." >&2
        exit 1
    fi
    # A worker killed by a segfault / memory exhaustion: exit 139 so the nf-core retry strategy re-runs the task
    # with more memory and half as many workers
    if grep -qE "TerminatedWorkerError|SIGSEGV|SIGKILL|MemoryError" ${prefix}.allo.log; then
        echo "ERROR: an Allo worker crashed (usually memory exhaustion with $workers workers on attempt ${task.attempt})." >&2
        echo "       Retrying with fewer workers / more memory; lower --allo_max_workers or raise memory if this persists." >&2
        exit 139
    fi
    if [ \$ALLO_RC -ne 0 ]; then
        echo "ERROR: Allo exited with status \$ALLO_RC, see ${prefix}.allo.log" >&2
        exit \$ALLO_RC
    fi

    # Allo exits 0 on some fatal errors; make sure it actually produced output
    if [ ! -s ${prefix}.allo.sam ] || ! grep -q "Allocation finished" ${prefix}.allo.log; then
        echo "ERROR: Allo did not finish successfully, see ${prefix}.allo.log" >&2
        exit 1
    fi

    # Compress Allo's SAM output to BAM with pysam (always present: Allo depends on it) and drop the SAM
    ALLO_PY=\$(head -n 1 "\$(command -v allo)" | sed 's/^#![[:space:]]*//')
    \$ALLO_PY -c "import sys, pysam; i = pysam.AlignmentFile(sys.argv[1]); o = pysam.AlignmentFile(sys.argv[2], 'wb', template=i); [o.write(r) for r in i]; o.close(); i.close()" ${prefix}.allo.sam ${prefix}.allo.bam
    rm -f ${prefix}.allo.sam

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
    touch ${prefix}.allo.bam ${prefix}.allo.log
    printf "${meta.id}\\treads\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}.allo_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allo: 1.2.0
    END_VERSIONS
    """
}
