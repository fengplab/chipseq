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
    path  sizes
    path  alias   // optional chromosome alias table ([] if none)

    output:
    tuple val(meta), path("*.features.bed")    , emit: bed
    tuple val(meta), path("*.chrom_report.tsv"), emit: report
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
                if (h == "thickstart") ts = i
                if (h == "thickend") te = i
                if (h == "blocksizes") bs = i
                if (h == "description") ds = i
            }
            # UCSC bigRmsk (bigRmskBed): chromStart/chromEnd are the *visualisation* span, which can include
            # unaligned parts of the repeat consensus; thickStart/thickEnd are the aligned genomic span
            rmsk = (ts && te && bs && ds) ? 1 : 0
            next
        }
        /^#/ { next }
        {
            # bigRmsk without a header line (older bigBedToBed): 14 fields and a 'name#class/family' name
            if (NR == 1 && !ts && NF == 14 && \$4 ~ /#/) { rmsk = 1; ts = 7; te = 8 }
            start = \$2; end = \$3
            if (rmsk && \$te ~ /^[0-9]+\$/ && \$ts ~ /^[0-9]+\$/ && \$te + 0 > \$ts + 0 && \$ts + 0 >= \$2 + 0 && \$te + 0 <= \$3 + 0) {
                start = \$ts; end = \$te; n_thick++
            }
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
            print \$1, start, end, name, score, strand, cls, fam
        }
        END { if (rmsk) printf "bigRmsk format detected: %d features use the aligned span (thickStart-thickEnd)\\n", n_thick > "/dev/stderr" }' raw.bed \\
        > normalised.bed

    # Harmonise chromosome names with the genome: exact match, then an optional alias table (any name on a line
    # maps to the name on that line that exists in the genome, e.g. UCSC chromAlias.txt), then adding/removing a
    # 'chr' prefix (plus chrM <-> MT). Features on chromosomes still absent from the genome are dropped and reported.
    awk -v alias_file="${alias}" -v report=${prefix}.chrom_report.tsv -v strict=${params.chrom_names_strict ? 1 : 0} 'BEGIN { FS = OFS = "\\t" }
        FILENAME == ARGV[1] { g[\$1] = 1; ng++; next }
        FILENAME == alias_file {
            if (\$0 ~ /^#/) next
            n = split(\$0, a, /[ \\t]+/); tgt = ""
            for (i = 1; i <= n; i++) if (a[i] in g) { tgt = a[i]; break }
            if (tgt != "") for (i = 1; i <= n; i++) if (a[i] != "") al[a[i]] = tgt
            next
        }
        {
            c = \$1; t = ""
            if (c in g)                                    { t = c; how = "exact" }
            else if (c in al)                              { t = al[c]; how = "alias" }
            else if (("chr" c) in g)                       { t = "chr" c; how = "add_chr" }
            else if (c ~ /^chr/ && (substr(c, 4) in g))    { t = substr(c, 4); how = "remove_chr" }
            else if ((c == "chrM" || c == "M") && ("MT" in g))   { t = "MT"; how = "mito" }
            else if (c == "MT" && ("chrM" in g))           { t = "chrM"; how = "mito" }
            if (t == "") { dropped[c]++; nd++; next }
            if (!((c SUBSEP t) in seen)) { seen[c SUBSEP t] = 1; pairs[++np] = c "\\t" t "\\t" how }
            kept[c]++; nk++
            \$1 = t; print
        }
        END {
            print "feature_chrom", "genome_chrom", "match", "n_features" > report
            for (i = 1; i <= np; i++) { split(pairs[i], p, "\\t"); print p[1], p[2], p[3], kept[p[1]] >> report }
            for (c in dropped) print c, "NA", "not_in_genome", dropped[c] >> report
            printf "Chromosome matching: %d features kept, %d dropped (not in genome)\\n", nk, nd > "/dev/stderr"
            if (nk == 0) {
                printf "%s: none of the feature chromosome names match the genome FASTA.\\n", (strict ? "ERROR" : "WARNING") > "/dev/stderr"
                k = 0; printf "  feature chromosomes, e.g.:" > "/dev/stderr"; for (c in dropped) { if (k++ < 5) printf " %s", c > "/dev/stderr" }
                k = 0; printf "\\n  genome chromosomes,  e.g.:" > "/dev/stderr"; for (c in g) { if (k++ < 5) printf " %s", c > "/dev/stderr" }
                printf "\\n  Provide a chromosome alias table with --feature_chrom_alias (e.g. UCSC <assembly>.chromAlias.txt).\\n" > "/dev/stderr"
                if (strict) exit 1
                printf "  This feature set will be skipped (set --chrom_names_strict true to stop the run instead).\\n" > "/dev/stderr"
            }
        }' $sizes ${alias ?: '/dev/null'} normalised.bed \\
        | LC_ALL=C sort -k1,1 -k2,2n > ${prefix}.features.bed

    rm -f normalised.bed

    rm -f raw.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
