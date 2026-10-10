/*
 * Check that a pre-built Bowtie2 index contains every sequence of the genome FASTA in use
 * (e.g. that an index given together with --viral_fasta was built from the combined host + viral genome)
 */
process BOWTIE2_INDEX_CHECK {
    tag "$fai"
    label 'process_single'

    conda "bioconda::bowtie2=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.2--py39h6fed5c7_0' :
        'biocontainers/bowtie2:2.5.2--py39h6fed5c7_0' }"

    input:
    tuple val(meta), path(index)
    path fai

    output:
    path "index_check.txt", emit: report
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found in $index" 1>&2 && exit 1

    bowtie2-inspect -n "\$INDEX" | awk '{print \$1}' | sort > index_names.txt
    cut -f1 $fai | sort > genome_names.txt
    comm -13 index_names.txt genome_names.txt > missing.txt
    {
        echo "index sequences : \$(wc -l < index_names.txt)"
        echo "genome sequences: \$(wc -l < genome_names.txt)"
        echo "missing in index: \$(wc -l < missing.txt)"
        cat missing.txt
    } > index_check.txt
    if [ -s missing.txt ]; then
        echo "ERROR: the pre-built --bowtie2_index does not contain \$(wc -l < missing.txt) sequence(s) of the genome in use:" 1>&2
        head -n 20 missing.txt 1>&2
        echo "Build the index from the combined host + viral FASTA (run once without --bowtie2_index and --save_reference, then reuse genome/index/bowtie2)." 1>&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
