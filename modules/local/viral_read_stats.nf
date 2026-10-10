/*
 * Count and extract reads on viral genomes (--viral_fasta) in a duplicate-marked BAM aligned to host + viruses.
 *
 * Per sample and virus: reads, reads after removing duplicates, % of all mapped reads, reads per million mapped,
 * mean depth and % of the viral genome covered. Counts are reads (both mates of a pair count), primary alignments
 * only, QC-fail removed. Viral reads are also written to <sample>.viral.bam for viewing (e.g. in IGV).
 */
process VIRAL_READ_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat)
    path contigs

    output:
    tuple val(meta), path("*.viral.bam"), path("*.viral.bam.bai"), emit: bam
    tuple val(meta), path("*.viral_reads.tsv")                     , emit: tsv
    tuple val(meta), path("*.viral_reads_mqc.tsv")                 , emit: mqc
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Viral sequences as BED (whole sequences)
    awk -F'\\t' 'NF >= 2 && \$1 != "" { print \$1 "\\t0\\t" \$2 }' $contigs > viral.bed

    # All mapped reads in the library (host + viral): primary, QC-passed, from the duplicate-marked flagstat
    mapped=\$(awk '/ primary mapped /    { print \$1; exit }' $flagstat)
    dups=\$(awk '/ primary duplicates / { print \$1; exit }' $flagstat)
    if [ -z "\$mapped" ] || [ -z "\$dups" ]; then
        mapped=\$(samtools view -c -@ $task.cpus -F 0xB04 $bam)
        dups=\$(samtools view -c -@ $task.cpus -F 0xB04 -f 0x400 $bam)
    fi
    mapped_dedup=\$(( mapped - dups ))

    # Reads on the viral sequences
    samtools view -@ $task.cpus -b -M -L viral.bed -o ${prefix}.viral.bam $bam
    samtools index ${prefix}.viral.bam

    # Depth and breadth (duplicates, secondary, QC-fail and unmapped excluded: samtools coverage defaults)
    samtools coverage ${prefix}.viral.bam > coverage.txt

    printf "sample\\tvirus\\tlength_bp\\treads\\treads_dedup\\tpct_of_mapped_reads\\treads_per_million_dedup\\tmean_depth_dedup\\tpct_genome_covered\\tdescription\\n" > ${prefix}.viral_reads.tsv
    : > ${prefix}.viral_reads_mqc.tsv
    while IFS=\$'\\t' read -r name length desc; do
        [ -n "\$name" ] || continue
        printf "%s\\t0\\t%s\\n" "\$name" "\$length" > one.bed
        n=\$(samtools view -c -F 0xB04 -M -L one.bed ${prefix}.viral.bam)
        nd=\$(samtools view -c -F 0xF04 -M -L one.bed ${prefix}.viral.bam)
        read -r depth covered < <(awk -F'\\t' -v c="\$name" '\$1 == c { print \$7, \$6 }' coverage.txt)
        awk -v s="${meta.id}" -v v="\$name" -v L="\$length" -v n="\$n" -v nd="\$nd" -v m="\$mapped" -v md="\$mapped_dedup" \\
            -v dp="\${depth:-0}" -v cv="\${covered:-0}" -v d="\$desc" -v mqc="${prefix}.viral_reads_mqc.tsv" 'BEGIN {
            pct = (m > 0) ? 100 * n / m : 0
            rpm = (md > 0) ? 1e6 * nd / md : 0
            printf "%s\\t%s\\t%d\\t%d\\t%d\\t%.4f\\t%.2f\\t%.3f\\t%.2f\\t%s\\n", s, v, L, n, nd, pct, rpm, dp, cv, d
            printf "%s | %s\\t%s\\t%s\\t%d\\t%d\\t%.4f\\t%.2f\\t%.3f\\t%.2f\\n", s, v, s, v, n, nd, pct, rpm, dp, cv >> mqc
        }' >> ${prefix}.viral_reads.tsv
    done < $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.viral.bam ${prefix}.viral.bam.bai ${prefix}.viral_reads.tsv ${prefix}.viral_reads_mqc.tsv
    echo '"${task.process}": {samtools: 1.20}' > versions.yml
    """
}
