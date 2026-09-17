include { SAMTOOLS_VIEW_MULTIMAPPED } from '../../modules/local/samtools_view_multimapped'
include { BEDTOOLS_COVERAGE_REPEATS } from '../../modules/local/bedtools_coverage_repeats'

//
// SUBWORKFLOW: Separate multi-mapped reads from a filtered/deduplicated BAM file and
//              quantify their coverage over annotated repetitive genome elements.
//
workflow BAM_MULTIMAP_REPEATS {
    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    ch_repeat_bed // channel: [ path(bed) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Split reads below the MAPQ threshold (multi-mapped) out of the BAM file
    //
    SAMTOOLS_VIEW_MULTIMAPPED (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MULTIMAPPED.out.versions.first())

    //
    // MODULE: Quantify multi-mapped read coverage over repetitive elements
    //
    BEDTOOLS_COVERAGE_REPEATS (
        SAMTOOLS_VIEW_MULTIMAPPED.out.bam,
        ch_repeat_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE_REPEATS.out.versions.first())

    emit:
    multimapped_bam = SAMTOOLS_VIEW_MULTIMAPPED.out.bam       // channel: [ val(meta), path(bam) ]
    unique_bam       = SAMTOOLS_VIEW_MULTIMAPPED.out.unique_bam // channel: [ val(meta), path(bam) ]
    coverage         = BEDTOOLS_COVERAGE_REPEATS.out.coverage   // channel: [ val(meta), path(txt) ]
    summary          = BEDTOOLS_COVERAGE_REPEATS.out.summary    // channel: [ val(meta), path(tsv) ]
    versions         = ch_versions                              // channel: [ versions.yml ]
}
