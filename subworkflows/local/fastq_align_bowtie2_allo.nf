//
// Alignment with Bowtie2 reporting multi-mapped reads (-k N), multi-mapped read allocation with Allo,
// then coordinate sorting, indexing and samtools stats on the Allo output.
//
// raw reads -> BOWTIE2_ALIGN (-k N, unsorted) -> SAMTOOLS_COLLATE (SAM, grouped by read name)
//           -> ALLO (--mixed) -> ALLO_SAM_TO_BAM (clear 0x100 on allocated reads) -> BAM_SORT_STATS_SAMTOOLS
//

include { BOWTIE2_ALIGN           } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_COLLATE        } from '../../modules/local/samtools_collate'
include { ALLO                    } from '../../modules/local/allo'
include { ALLO_SAM_TO_BAM         } from '../../modules/local/allo_sam_to_bam'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BOWTIE2_ALLO {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_index          // channel: /path/to/bowtie2/index/
    save_unaligned    // val
    ch_fasta          // channel: [ val(meta), fasta ]
    allo_header       // file : MultiQC custom-content header for the Allo table

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2 keeping up to N alignments per read. Output is left unsorted
    // (Bowtie2 emits all alignments of a read/pair contiguously).
    //
    BOWTIE2_ALIGN ( ch_reads, ch_index, ch_fasta, save_unaligned, false )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Group by read name and write the raw alignments as SAM (Allo checks for a 'samtools collate' @PG line)
    //
    SAMTOOLS_COLLATE ( BOWTIE2_ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions.first())

    //
    // Allocate multi-mapped reads with Allo
    //
    ALLO ( SAMTOOLS_COLLATE.out.sam )
    ch_versions = ch_versions.mix(ALLO.out.versions.first())

    //
    // Make allocated reads primary and convert to BAM
    //
    ALLO_SAM_TO_BAM ( ALLO.out.sam )
    ch_versions = ch_versions.mix(ALLO_SAM_TO_BAM.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( ALLO_SAM_TO_BAM.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // Collate per-library Allo summaries into a single MultiQC custom-content table
    //
    ALLO
        .out
        .mqc
        .map { meta, tsv -> tsv }
        .collectFile(name: 'allo_allocation_mqc.tsv', seed: allo_header.text, sort: true, storeDir: "${params.outdir}/${params.aligner}/library/allo")
        .set { ch_allo_multiqc }

    emit:
    bam_orig         = BOWTIE2_ALIGN.out.bam                 // channel: [ val(meta), bam ]
    log_out          = BOWTIE2_ALIGN.out.log                 // channel: [ val(meta), log ]
    fastq            = BOWTIE2_ALIGN.out.fastq               // channel: [ val(meta), fastq ]
    allo_sam         = ALLO.out.sam                          // channel: [ val(meta), sam ]
    allo_log         = ALLO.out.log                          // channel: [ val(meta), log ]
    allo_multiqc     = ch_allo_multiqc                       // channel: path(allo_allocation_mqc.tsv)

    bam              = BAM_SORT_STATS_SAMTOOLS.out.bam       // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_STATS_SAMTOOLS.out.bai       // channel: [ val(meta), [ bai ] ]
    csi              = BAM_SORT_STATS_SAMTOOLS.out.csi       // channel: [ val(meta), [ csi ] ]
    stats            = BAM_SORT_STATS_SAMTOOLS.out.stats     // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats  // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                           // channel: [ versions.yml ]
}
