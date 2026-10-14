/*
 * Prepare genome intervals for filtering by removing regions in blacklist file.
 * Whole sequences listed in 'exclude' (viral genomes from --viral_fasta, names from --exclude_contigs) are left out
 * entirely, so the filtered BAM - and everything downstream of it - contains no reads on them.
 */
process GENOME_BLACKLIST_REGIONS {
    tag "$sizes"

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    path sizes
    path blacklist
    path exclude

    output:
    path '*.include_regions.bed', emit: bed
    path 'excluded_contigs.txt' , emit: excluded, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    def drop = exclude ? """
        # Leave out whole sequences (viral genomes, --exclude_contigs)
        cut -f1 $exclude | sed '/^[[:space:]]*\$/d' | sort -u > exclude_names.txt
        awk 'NR==FNR { x[\$1]=1; next } !(\$1 in x)' exclude_names.txt $sizes > sizes.kept
        awk 'NR==FNR { x[\$1]=1; next } (\$1 in x) { print \$1 "\\t" \$2 }' exclude_names.txt $sizes > excluded_contigs.txt
        missing=\$(awk 'NR==FNR { g[\$1]=1; next } !(\$1 in g) { m = m (m == "" ? "" : ",") \$1 } END { print m }' $sizes exclude_names.txt)
        if [ -n "\$missing" ]; then
            echo "ERROR: sequence(s) to exclude not found in the genome (check --exclude_contigs against the FASTA names): \$missing" 1>&2
            exit 1
        fi
        [ -s sizes.kept ] || { echo "ERROR: every genome sequence was excluded" 1>&2; exit 1; }
        """ : "cp $sizes sizes.kept"
    if (blacklist) {
        """
        $drop
        # Blacklist intervals on excluded sequences are dropped too
        awk 'NR==FNR { k[\$1]=1; next } (\$1 in k)' sizes.kept $blacklist > blacklist.kept
        sortBed -i blacklist.kept -g sizes.kept | complementBed -i stdin -g sizes.kept > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        $drop
        awk '{print \$1, '0' , \$2}' OFS='\t' sizes.kept > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }

    stub:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    """
    touch $file_out
    ${exclude ? 'touch excluded_contigs.txt' : ''}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.30.0
    END_VERSIONS
    """
}
