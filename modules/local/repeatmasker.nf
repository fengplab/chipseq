/*
 * Run RepeatMasker on the reference genome (only used with --run_repeatmasker when no pre-generated
 * RepeatMasker BigBed/annotation is supplied). This is slow for large genomes: publish the resulting
 * BigBed (genome/annotation/repeatmasker.bb) and pass it via --repeatmasker_bigbed on subsequent runs.
 */
process REPEATMASKER {
    tag "$fasta"
    label 'process_high'
    label 'process_long'

    conda "bioconda::repeatmasker=4.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    path fasta
    path lib

    output:
    path "*.rmsk.out"   , emit: out
    path "*.rmsk.tbl"   , emit: tbl, optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def lib_arg = lib ? "-lib $lib" : ''
    def prefix  = fasta.getBaseName(fasta.name.endsWith('.gz') ? 2 : 1)
    def unzip   = fasta.name.endsWith('.gz') ? "gunzip -c $fasta > genome.fa" : "ln -s $fasta genome.fa"
    """
    $unzip

    RepeatMasker \\
        $lib_arg \\
        -pa $task.cpus \\
        -dir rm_out \\
        $args \\
        genome.fa

    mv rm_out/genome.fa.out ${prefix}.rmsk.out
    [ -f rm_out/genome.fa.tbl ] && mv rm_out/genome.fa.tbl ${prefix}.rmsk.tbl || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v 2>&1 | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
