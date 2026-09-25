/*
 * Convert a (pre-generated or freshly built) BigBed of features to a normalised, sorted BED8:
 *   chrom  start  end  name  score  strand  class  family
 * 'class'/'family' are taken from autoSql fields named repClass/class/type and repFamily/family when present,
 * from RepeatMasker style names ('AluY#SINE/Alu'), or otherwise derived from the feature name
 * (CenSat: 'hsat2_2(...)' -> 'hsat2').
 */
process BIGBED_TO_BED {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ucsc-bigbedtobed=482"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bigbedtobed:482--h0b57e2e_0' :
        'biocontainers/ucsc-bigbedtobed:482--h0b57e2e_0' }"

    input:
    tuple val(meta), path(bigbed)

    output:
    tuple val(meta), path("*.features.bed"), emit: bed
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '482' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bigBedToBed -header $bigbed raw.bed 2>/dev/null || bigBedToBed $bigbed raw.bed

    awk -v type="${meta.id}" 'BEGIN { FS = OFS = "\\t"; cc = 0; fc = 0 }
        NR == 1 && /^#/ {
            sub(/^#/, "")
            for (i = 1; i <= NF; i++) {
                h = tolower(\$i)
                if (!cc && (h == "repclass" || h == "class" || h == "type" || h == "classname")) cc = i
                if (!fc && (h == "repfamily" || h == "family")) fc = i
            }
            next
        }
        /^#/ { next }
        {
            name   = (NF >= 4 && \$4 != "") ? \$4 : \$1 ":" \$2 "-" \$3
            score  = (NF >= 5 && \$5 ~ /^[0-9.]+\$/) ? int(\$5) : 0
            strand = (NF >= 6 && (\$6 == "+" || \$6 == "-")) ? \$6 : "."
            cls = ""; fam = ""
            if (cc && cc <= NF) cls = \$cc
            if (fc && fc <= NF) fam = \$fc
            if (cls == "" && index(name, "#") > 0) {
                split(name, p, "#"); name = p[1]; cf = p[2]
                n = split(cf, q, "/"); cls = q[1]; fam = (n > 1) ? q[2] : q[1]
            }
            if (cls == "") {
                cls = name
                if (type == "censat") { sub(/\\(.*\$/, "", cls); while (cls ~ /_[0-9]+\$/) sub(/_[0-9]+\$/, "", cls) }
            }
            if (cls == "") cls = "NA"
            if (fam == "") fam = cls
            gsub(/ /, "_", cls); gsub(/ /, "_", fam)
            print \$1, \$2, \$3, name, score, strand, cls, fam
        }' raw.bed \\
        | LC_ALL=C sort -k1,1 -k2,2n > ${prefix}.features.bed

    rm -f raw.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
