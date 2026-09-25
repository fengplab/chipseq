//
// Motif discovery on called peaks (per-sample and consensus) with XSTREME.
// --maxw is set to the length of the shortest read in the final filtered, sorted alignment file(s)
// (the *.mLb.clN.sorted.bam files that were used for peak calling), taken from the samtools stats
// 'RL' (read length) histogram of those files.
//

include { PEAKS_TO_FASTA } from '../../modules/local/peaks_to_fasta'
include { MEME_XSTREME   } from '../../modules/local/meme_xstreme'

//
// Shortest read length from a samtools stats file ('RL <length> <count>' lines)
//
def getMinReadLength(stats_file) {
    def min_len = null
    stats_file.eachLine { line ->
        if (line.startsWith('RL\t')) {
            def fields = line.split('\t')
            def len    = fields[1].toInteger()
            def count  = fields[2].toLong()
            if (count > 0 && (min_len == null || len < min_len)) {
                min_len = len
            }
        }
    }
    return min_len
}

//
// Apply the optional STREME-compatible cap to the requested motif width
//
def resolveMaxw(id, min_len, cap) {
    if (min_len == null) {
        log.warn "[XSTREME] Could not determine the shortest read length for '${id}'; falling back to the XSTREME default --maxw (15)."
        return 15
    }
    if (cap && cap > 0 && min_len > cap) {
        log.warn "[XSTREME] Shortest read for '${id}' is ${min_len} bp but STREME only supports motif widths <= ${cap}; using --maxw ${cap}. Set '--xstreme_maxw_cap 0' to pass ${min_len} unchanged."
        return cap
    }
    return min_len
}

workflow PEAKS_MOTIFS_XSTREME {
    take:
    ch_peaks            // channel: [ val(meta), peaks ]             per-sample MACS3 peaks (meta.antibody set)
    ch_consensus_peaks  // channel: [ val(meta), bed ]               consensus peaks, meta.id == antibody (may be empty)
    ch_bam_stats        // channel: [ val(meta), stats ]             samtools stats of the final *.mLb.clN.sorted.bam
    ch_fasta            // channel: path(fasta)
    ch_fai              // channel: path(fai)
    ch_motif_db         //    path: optional known motif database(s) in MEME format ([] if none)
    maxw_cap            // integer: cap applied to --maxw (0 = no cap)

    main:

    ch_versions = Channel.empty()

    // [ sample_id, min_read_length ]
    ch_bam_stats
        .map { meta, stats -> [ meta.id, getMinReadLength(stats) ] }
        .set { ch_min_len }

    // Per-sample peaks: [ meta, peaks, maxw ]
    ch_peaks
        .map { meta, peaks -> [ meta.id, meta, peaks ] }
        .join(ch_min_len)
        .map { id, meta, peaks, min_len -> [ meta, peaks, resolveMaxw(id, min_len, maxw_cap) ] }
        .set { ch_sample_peaks_maxw }

    // Consensus peaks: shortest read across all IP samples of that antibody
    ch_peaks
        .map { meta, peaks -> [ meta.id, meta.antibody ] }
        .join(ch_min_len)
        .map { id, antibody, min_len -> [ antibody, min_len ] }
        .groupTuple()
        .map { antibody, lens -> [ antibody, lens.findAll { it != null }.min() ] }
        .set { ch_antibody_min_len }

    ch_consensus_peaks
        .map { meta, bed -> [ meta.id, meta, bed ] }
        .join(ch_antibody_min_len)
        .map { antibody, meta, bed, min_len ->
            def meta_new = meta + [ id: "${antibody}.consensus_peaks".toString(), antibody: antibody, consensus: true ]
            [ meta_new, bed, resolveMaxw(antibody, min_len, maxw_cap) ]
        }
        .set { ch_consensus_peaks_maxw }

    //
    // Extract peak sequences
    //
    PEAKS_TO_FASTA (
        ch_sample_peaks_maxw.mix(ch_consensus_peaks_maxw),
        ch_fasta.collect(),
        ch_fai.collect()
    )
    ch_versions = ch_versions.mix(PEAKS_TO_FASTA.out.versions.first())

    //
    // Run XSTREME
    //
    MEME_XSTREME (
        PEAKS_TO_FASTA.out.fasta,
        ch_motif_db
    )
    ch_versions = ch_versions.mix(MEME_XSTREME.out.versions.first())

    emit:
    fasta    = PEAKS_TO_FASTA.out.fasta     // channel: [ val(meta), fasta, maxw ]
    results  = MEME_XSTREME.out.results     // channel: [ val(meta), dir ]
    html     = MEME_XSTREME.out.html        // channel: [ val(meta), html ]
    motifs   = MEME_XSTREME.out.motifs      // channel: [ val(meta), meme ]
    versions = ch_versions                  // channel: [ versions.yml ]
}
