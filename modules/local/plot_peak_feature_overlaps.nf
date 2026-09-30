/*
 * Per-sample plots + MultiQC table of peaks covering RepeatMasker / CenSat features
 */
process PLOT_PEAK_FEATURE_OVERLAPS {
    tag "$feature.id"
    label 'process_low'

    conda "conda-forge::r-base conda-forge::r-optparse conda-forge::r-ggplot2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(feature), path(summaries)
    val   criterion   // human-readable overlap filter description

    output:
    path "*.class_summary.combined.tsv", emit: tsv
    path "*.pdf"                       , emit: pdf
    path "*_mqc.tsv"                   , emit: multiqc
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    """
    plot_peak_feature_overlaps.r \\
        --summary_files ${summaries.collect{ it.toString() }.sort().join(',')} \\
        --feature_set ${feature.id} \\
        --criterion '${criterion}' \\
        --outprefix ${feature.id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
