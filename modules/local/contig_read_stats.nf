/*
 * Count reads on the whole sequences that are left out of the filtered BAM: viral genomes (--viral_fasta) and
 * --exclude_contigs (e.g. chrM). Counted in the duplicate-marked merged-library BAM, before those sequences are removed.
 *
 * Per sample and sequence: reads, reads after removing duplicates, and the de-duplicated reads as a percentage of
 * all reads given to the aligner after removing duplicates:
 *
 *     pct_of_total_reads_dedup = 100 * reads_dedup / (total_reads - duplicates)
 *
 *     total_reads : reads given to Bowtie2 (after trimming), summed over every sequencing run merged into this library,
 *                   read from the Bowtie2 logs (2 x pairs + unpaired reads). Taken from the logs because Allo drops
 *                   unaligned reads, so the BAM itself no longer contains them.
 *     duplicates  : primary duplicate reads in the duplicate-marked BAM (samtools flagstat). Duplicates are only marked
 *                   among aligned reads, so this is the number of reads removed by de-duplication.
 *
 * Reads are counted as single reads (both mates of a pair count), primary alignments only, QC-failed and
 * supplementary excluded. With --with_allo, multi-mapped reads count where Allo placed them.
 * Reads on viral genomes are also written to <sample>.viral.bam (e.g. for IGV). No BAM is written for other sequences.
 */
process CONTIG_READ_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat), path(bowtie2_logs)
    path excluded   // excluded_contigs.txt from GENOME_BLACKLIST_REGIONS: name, length (viral + --exclude_contigs)
    path viral      // viral_contigs.tsv from COMBINE_VIRAL_GENOME: name, length, description ([] without --viral_fasta)

    output:
    tuple val(meta), path("*.viral.bam"), path("*.viral.bam.bai"), emit: bam, optional: true
    tuple val(meta), path("*.contig_reads.tsv")                    , emit: tsv
    tuple val(meta), path("*.contig_reads_mqc.tsv")                , emit: mqc
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def get_viral = viral ? "cp $viral viral.tsv" : ": > viral.tsv"
    """
    set -euo pipefail
    fail() { echo "ERROR [${meta.id}]: \$*" >&2; exit 1; }

    # ---- Sequences to report -------------------------------------------------------------------------------
    # contigs.tsv: name, length, type ('viral' if in viral_contigs.tsv, otherwise 'excluded'), description
    $get_viral
    awk -F'\\t' -v OFS='\\t' '
        FILENAME == "viral.tsv" { if (\$1 != "") { v[\$1] = 1; d[\$1] = \$3 }; next }
        NF >= 2 && \$1 != ""     { print \$1, \$2, ((\$1 in v) ? "viral" : "excluded"), ((\$1 in d) ? d[\$1] : "") }
    ' viral.tsv $excluded > contigs.tsv
    [ -s contigs.tsv ] || fail "no sequences to report in $excluded"
    n_contigs=\$(wc -l < contigs.tsv)

    # Every sequence must be in the BAM header under the same name
    samtools view -H $bam | awk -F'\\t' '\$1 == "@SQ" { for (i = 2; i <= NF; i++) if (\$i ~ /^SN:/) print substr(\$i, 4) }' > bam_contigs.txt
    missing=\$(awk 'FILENAME == "bam_contigs.txt" { h[\$1] = 1; next } !(\$1 in h) { m = m (m == "" ? "" : ",") \$1 } END { print m }' bam_contigs.txt contigs.tsv)
    [ -z "\$missing" ] || fail "sequence(s) not in the BAM header: \$missing"

    # BED of whole sequences. BED is 0-based half-open, so [0, length) covers every base of a sequence whose length
    # comes from samtools faidx (a base count). samtools view -L uses the BED as given; no conversion is needed.
    awk -F'\\t' -v OFS='\\t' '{ print \$1, 0, \$2 }' contigs.tsv > contigs.bed
    awk -F'\\t' -v OFS='\\t' '\$3 == "viral" { print \$1, 0, \$2 }' contigs.tsv > viral.bed

    # ---- Library totals --------------------------------------------------------------------------------------
    # Reads given to Bowtie2, from each run's log. First line: '<N> reads; of these:' (N = pairs + unpaired reads),
    # then '  <P> (x%) were paired; of these:' and/or '  <U> (x%) were unpaired; of these:'. Reads = 2P + U.
    total=0
    n_logs=0
    for log in $bowtie2_logs; do
        n=\$(awk '/^[0-9]+ reads; of these:\$/ { print \$1; exit }' "\$log")
        p=\$(awk '/^  [0-9]+ \\([0-9.]+%\\) were paired; of these:\$/ { print \$1; exit }' "\$log")
        u=\$(awk '/^  [0-9]+ \\([0-9.]+%\\) were unpaired; of these:\$/ { print \$1; exit }' "\$log")
        [ -n "\$n" ] || fail "no '<N> reads; of these:' line in Bowtie2 log \$log"
        p=\${p:-0}; u=\${u:-0}
        [ \$(( p + u )) -eq "\$n" ] || fail "Bowtie2 log \$log: paired (\$p) + unpaired (\$u) != reads (\$n)"
        reads=\$(( 2 * p + u ))
        echo "[contig_read_stats] ${meta.id}: Bowtie2 log \$log: \$p pairs + \$u unpaired = \$reads reads" >&2
        total=\$(( total + reads ))
        n_logs=\$(( n_logs + 1 ))
    done
    [ "\$n_logs" -gt 0 ] || fail "no Bowtie2 log for this library"
    [ "\$total" -gt 0 ] || fail "Bowtie2 logs report 0 reads"

    # Primary, QC-passed mapped reads and duplicates in the duplicate-marked BAM (first column of samtools flagstat)
    mapped=\$(awk '/ primary mapped \\(/  { print \$1; exit }' $flagstat)
    dups=\$(awk '/ primary duplicates\$/ { print \$1; exit }' $flagstat)
    if [ -z "\$mapped" ] || [ -z "\$dups" ]; then
        mapped=\$(samtools view -c -@ $task.cpus -F 0xB04 $bam)
        dups=\$(samtools view -c -@ $task.cpus -F 0xB04 -f 0x400 $bam)
    fi
    [ "\$mapped" -le "\$total" ] || fail "mapped reads in the BAM (\$mapped) exceed reads given to Bowtie2 (\$total): the Bowtie2 logs do not match this library"
    [ "\$dups" -le "\$mapped" ]  || fail "duplicates (\$dups) exceed mapped reads (\$mapped)"
    total_dedup=\$(( total - dups ))
    echo "[contig_read_stats] ${meta.id}: \$n_logs Bowtie2 log(s), \$total reads in, \$dups duplicates, \$total_dedup reads after de-duplication, \$mapped mapped (incl. duplicates)" >&2

    # ---- Per-sequence counts ----------------------------------------------------------------------------------
    # Reads on the reported sequences only (temporary; deleted below)
    samtools view -@ $task.cpus -b -M -L contigs.bed -o targets.bam $bam
    samtools index targets.bam

    # Depth and breadth: samtools coverage defaults skip unmapped, secondary, QC-failed and duplicate reads.
    # -d 0 removes the default depth cap (1,000,000), which a high-copy sequence could reach.
    samtools coverage -d 0 targets.bam > coverage.txt

    printf "sample\\tcontig\\tcontig_type\\tlength_bp\\ttotal_reads\\tduplicates\\ttotal_reads_dedup\\tmapped_reads\\treads\\treads_dedup\\tpct_of_total_reads_dedup\\tmean_depth_dedup\\tpct_genome_covered\\tdescription\\n" > ${prefix}.contig_reads.tsv
    : > ${prefix}.contig_reads_mqc.tsv
    sum_reads=0
    while IFS=\$'\\t' read -r name length type desc; do
        printf "%s\\t0\\t%s\\n" "\$name" "\$length" > one.bed
        n=\$(samtools view -c -F 0xB04 -M -L one.bed targets.bam)     # primary, mapped, not QC-failed/supplementary
        nd=\$(samtools view -c -F 0xF04 -M -L one.bed targets.bam)    # ... and not duplicates
        [ "\$nd" -le "\$n" ] || fail "\$name: de-duplicated reads (\$nd) exceed reads (\$n)"
        sum_reads=\$(( sum_reads + n ))
        depth=\$(awk -F'\\t' -v c="\$name" '\$1 == c { print \$7 }' coverage.txt)
        covered=\$(awk -F'\\t' -v c="\$name" '\$1 == c { print \$6 }' coverage.txt)
        [ -n "\$depth" ] && [ -n "\$covered" ] || fail "\$name missing from samtools coverage output"
        echo "[contig_read_stats] ${meta.id}: \$name (\$type): \$n reads, \$nd after de-duplication" >&2
        awk -v s="${meta.id}" -v c="\$name" -v t="\$type" -v L="\$length" -v tot="\$total" -v du="\$dups" -v td="\$total_dedup" \\
            -v m="\$mapped" -v n="\$n" -v nd="\$nd" -v dp="\$depth" -v cv="\$covered" -v d="\$desc" \\
            -v mqc="${prefix}.contig_reads_mqc.tsv" 'BEGIN {
            pct = 100 * nd / td
            printf "%s\\t%s\\t%s\\t%.0f\\t%.0f\\t%.0f\\t%.0f\\t%.0f\\t%.0f\\t%.0f\\t%.6f\\t%.3f\\t%.2f\\t%s\\n", s, c, t, L, tot, du, td, m, n, nd, pct, dp, cv, d
            printf "%s | %s\\t%s\\t%s\\t%s\\t%.0f\\t%.0f\\t%.6f\\t%.3f\\t%.2f\\n", s, c, s, c, t, n, nd, pct, dp, cv >> mqc
        }' >> ${prefix}.contig_reads.tsv
    done < contigs.tsv
    [ "\$sum_reads" -le "\$mapped" ] || fail "reads on the reported sequences (\$sum_reads) exceed mapped reads (\$mapped)"
    rows=\$(( \$(wc -l < ${prefix}.contig_reads.tsv) - 1 ))
    [ "\$rows" -eq "\$n_contigs" ] || fail "wrote \$rows rows for \$n_contigs sequences"

    # Viral reads as a BAM for viewing
    if [ -s viral.bed ]; then
        samtools view -@ $task.cpus -b -M -L viral.bed -o ${prefix}.viral.bam targets.bam
        samtools index ${prefix}.viral.bam
    fi
    rm -f targets.bam targets.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample\\tcontig\\tcontig_type\\tlength_bp\\ttotal_reads\\tduplicates\\ttotal_reads_dedup\\tmapped_reads\\treads\\treads_dedup\\tpct_of_total_reads_dedup\\tmean_depth_dedup\\tpct_genome_covered\\tdescription\\n" > ${prefix}.contig_reads.tsv
    touch ${prefix}.contig_reads_mqc.tsv
    if [ -n "$viral" ]; then touch ${prefix}.viral.bam ${prefix}.viral.bam.bai; fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
