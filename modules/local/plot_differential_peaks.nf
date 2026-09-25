/*
 * DESeq2 differential binding between sample groups on consensus peaks, and visualisation of the
 * differential peaks together with their HOMER and RepeatMasker/CenSat annotations
 */
process PLOT_DIFFERENTIAL_PEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(meta), path(counts), path(homer_annotation), val(feature_names), path(feature_files)

    output:
    tuple val(meta), path("*.results.txt")          , emit: results     , optional: true
    tuple val(meta), path("*.significant.txt")      , emit: significant , optional: true
    tuple val(meta), path("*.significant.bed")      , emit: bed         , optional: true
    tuple val(meta), path("*_enrichment.txt")       , emit: enrichment  , optional: true
    tuple val(meta), path("*.pdf")                  , emit: pdf         , optional: true
    tuple val(meta), path("*.differential_summary.txt"), emit: summary  , optional: true
    tuple val(meta), path("*.normalised_counts.txt"), emit: norm_counts , optional: true
    tuple val(meta), path("*.rds")                  , emit: rds         , optional: true
    path  "*_mqc.tsv"                               , emit: multiqc     , optional: true
    path  "*.log"                                   , emit: log
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def homer  = homer_annotation ? "--homer_file ${homer_annotation}" : ''
    def ffiles = feature_files ? [feature_files].flatten() : []
    def fnames = feature_names ? [feature_names].flatten() : []
    def feats  = ffiles ? "--feature_files ${ffiles.join(',')} --feature_names ${fnames.join(',')}" : ''
    """
    plot_differential_peaks.r \\
        --count_file $counts \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $homer \\
        $feats \\
        $args \\
        2>&1 | tee ${prefix}.differential.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
