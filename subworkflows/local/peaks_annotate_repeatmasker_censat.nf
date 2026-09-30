//
// Annotate peaks with RepeatMasker and CenSat features.
//
// Feature sources, in order of precedence (per feature set):
//   RepeatMasker: --repeatmasker_bigbed (pre-generated, reused as-is)
//                 --repeatmasker_annotation (RepeatMasker .out, UCSC rmsk.txt[.gz] or BED -> BigBed)
//                 --run_repeatmasker (RepeatMasker is run on the genome -> BigBed)
//   CenSat      : --censat_bigbed (pre-generated, reused as-is)
//                 --censat_bed (CenSat BED -> BigBed)
// BigBeds that are built by the pipeline are published to <outdir>/genome/annotation/ so they can be passed
// back in with --repeatmasker_bigbed / --censat_bigbed on later runs and do not have to be re-calculated.
//

include { REPEATMASKER               } from '../../modules/local/repeatmasker'
include { FEATURES_TO_BED            } from '../../modules/local/features_to_bed'
include { UCSC_BEDTOBIGBED           } from '../../modules/local/ucsc_bedtobigbed'
include { BIGBED_TO_BED              } from '../../modules/local/bigbed_to_bed'
include { PEAK_FEATURE_INTERSECT     } from '../../modules/local/peak_feature_intersect'
include { PLOT_PEAK_FEATURE_OVERLAPS } from '../../modules/local/plot_peak_feature_overlaps'

workflow PEAKS_ANNOTATE_REPEATMASKER_CENSAT {
    take:
    ch_peaks           // channel: [ val(meta), peaks ]  per-sample MACS3 peaks
    ch_consensus_peaks // channel: [ val(meta), bed ]    consensus peaks (meta.id == antibody), may be empty
    ch_fasta           // channel: path(fasta)
    ch_chrom_sizes     // channel: path(chrom.sizes)
    min_overlap        // float  : minimum fraction of a feature that must be covered by a peak

    main:

    ch_versions   = Channel.empty()
    ch_bigbed     = Channel.empty() // [ [id:feature_set], bigbed ]
    ch_source     = Channel.empty() // [ [id:feature_set], annotation to convert ]
    ch_fasta_val  = ch_fasta.collect().map { it[0] }
    ch_sizes_val  = ch_chrom_sizes.collect().map { it[0] }
    ch_alias      = params.feature_chrom_alias ? file(params.feature_chrom_alias, checkIfExists: true) : []

    //
    // RepeatMasker
    //
    if (params.repeatmasker_bigbed) {
        ch_bigbed = ch_bigbed.mix(Channel.of([ [ id:'repeatmasker' ], file(params.repeatmasker_bigbed, checkIfExists: true) ]))
    } else if (params.repeatmasker_annotation) {
        ch_source = ch_source.mix(Channel.of([ [ id:'repeatmasker' ], file(params.repeatmasker_annotation, checkIfExists: true) ]))
    } else if (params.run_repeatmasker) {
        REPEATMASKER (
            ch_fasta_val,
            params.repeatmasker_lib ? file(params.repeatmasker_lib, checkIfExists: true) : []
        )
        ch_source   = ch_source.mix(REPEATMASKER.out.out.map { [ [ id:'repeatmasker' ], it ] })
        ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
    } else {
        log.warn "[Feature annotation] No RepeatMasker source given (--repeatmasker_bigbed, --repeatmasker_annotation or --run_repeatmasker); skipping RepeatMasker annotation."
    }

    //
    // CenSat
    //
    if (params.censat_bigbed) {
        ch_bigbed = ch_bigbed.mix(Channel.of([ [ id:'censat' ], file(params.censat_bigbed, checkIfExists: true) ]))
    } else if (params.censat_bed) {
        ch_source = ch_source.mix(Channel.of([ [ id:'censat' ], file(params.censat_bed, checkIfExists: true) ]))
    } else {
        log.warn "[Feature annotation] No CenSat source given (--censat_bigbed or --censat_bed); skipping CenSat annotation."
    }

    //
    // Build (and publish) BigBeds for sources that are not BigBed yet
    //
    FEATURES_TO_BED (
        ch_source,
        ch_sizes_val,
        ch_alias
    )
    ch_versions = ch_versions.mix(FEATURES_TO_BED.out.versions.first())

    UCSC_BEDTOBIGBED (
        FEATURES_TO_BED.out.bed,
        ch_sizes_val,
        file("$projectDir/assets/feature_bed6plus2.as", checkIfExists: true)
    )
    ch_bigbed   = ch_bigbed.mix(UCSC_BEDTOBIGBED.out.bigbed)
    ch_versions = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions.first())

    //
    // BigBed -> normalised BED8 used for the intersections
    //
    BIGBED_TO_BED (
        ch_bigbed,
        ch_sizes_val,
        ch_alias
    )
    ch_versions = ch_versions.mix(BIGBED_TO_BED.out.versions.first())

    //
    // Intersect every peak set (per-sample and consensus) with every feature set
    //
    ch_consensus_peaks
        .map {
            meta, bed ->
                [ meta + [ id: "${meta.id}.consensus_peaks".toString(), antibody: meta.id, consensus: true ], bed ]
        }
        .set { ch_consensus_meta }

    PEAK_FEATURE_INTERSECT (
        ch_peaks
            .map { meta, peaks -> [ meta + [ consensus: false ], peaks ] }
            .mix(ch_consensus_meta)
            .combine(BIGBED_TO_BED.out.bed),
        min_overlap
    )
    ch_versions = ch_versions.mix(PEAK_FEATURE_INTERSECT.out.versions.first())

    //
    // Per-sample summary plots + MultiQC tables (one per feature set)
    //
    PLOT_PEAK_FEATURE_OVERLAPS (
        PEAK_FEATURE_INTERSECT
            .out
            .summary
            .filter { meta, feature, tsv -> !meta.consensus }
            .map { meta, feature, tsv -> [ feature, tsv ] }
            .groupTuple(),
        min_overlap
    )
    ch_versions = ch_versions.mix(PLOT_PEAK_FEATURE_OVERLAPS.out.versions.first())

    emit:
    bigbed          = ch_bigbed                                   // channel: [ val(meta), bigbed ]
    features_bed    = BIGBED_TO_BED.out.bed                       // channel: [ val(meta), bed ]
    chrom_report    = BIGBED_TO_BED.out.report                    // channel: [ val(meta), tsv ]
    all             = PEAK_FEATURE_INTERSECT.out.all              // channel: [ val(meta), val(feature), tsv ]
    covered         = PEAK_FEATURE_INTERSECT.out.covered          // channel: [ val(meta), val(feature), tsv ]
    peak_annotation = PEAK_FEATURE_INTERSECT.out.peak_annotation  // channel: [ val(meta), val(feature), tsv ]
    summary         = PEAK_FEATURE_INTERSECT.out.summary          // channel: [ val(meta), val(feature), tsv ]
    plots           = PLOT_PEAK_FEATURE_OVERLAPS.out.pdf          // channel: [ pdf ]
    multiqc         = PLOT_PEAK_FEATURE_OVERLAPS.out.multiqc      // channel: [ tsv ]
    versions        = ch_versions                                 // channel: [ versions.yml ]
}
