#!/usr/bin/env python3
#SBATCH --job-name=preflight_reference
#SBATCH --partition defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/preflight_reference_%j.out
#SBATCH --error=logs/preflight_reference_%j.out
"""
Pre-flight check of the reference set-up before running nf-core/chipseq with the GenBank T2T assembly and
--viral_fasta. Checks, without running the pipeline, everything the reference steps will rely on:

  1. pipeline checkout   the local chipseq/ copy has the viral-genome version; params-file keys are known
  2. host FASTA          readable, unique names, contains the mitochondrion; not the RefSeq file without chrM
  3. --chrom_alias       every host sequence gets a name; no two names collide; mitochondrion -> chrM
  4. same assembly       (optional) chromosome lengths identical to the RefSeq FASTA the GTF was made for
  5. GTF                 chromosome names match the renamed genome; all genes lie inside the chromosomes
  6. --viral_fasta       readable, unique names, no clash with host names, sensible lengths and bases
  7. --exclude_contigs   every name exists in the final genome
  8. blacklist           chromosome names match; intervals inside the chromosomes
  9. --bowtie2_index     compatible with --chrom_alias / --viral_fasta (contents checked if bowtie2-inspect works)
 10. environment         Singularity images cached (prefetch_singularity_images.sh --check), free disk space

Reads the reference paths from nextflow-params.yml and the human GENOME/GTF lines of run_chipseq_slurm.sh;
override them with the options below. Reading the FASTA files takes a few minutes, so submit it as a job from
the launch directory:
    mkdir -p logs && sbatch preflight_reference.py
    python3 preflight_reference.py --fasta GCA_009914755.4_T2T-CHM13v2.0_genomic.fna      # or run directly
Exit status 0 = no FAIL lines.
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

RESULTS = {"PASS": 0, "WARN": 0, "FAIL": 0}
MITO_LEN = 16569


def say(level, msg):
    if level in RESULTS:
        RESULTS[level] += 1
    print(f"{level:5} {msg}", flush=True)


def section(title):
    print(f"\n== {title} ==", flush=True)


def open_any(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "r")


def short(names, n=8):
    names = list(names)
    return ", ".join(names[:n]) + (f" ... (+{len(names) - n} more)" if len(names) > n else "")


# ---------------------------------------------------------------------------------------------------------------
def read_params(path):
    """Flat 'key: value' YAML (as used for -params-file); comments and quotes stripped."""
    params = OrderedDict()
    if not path or not os.path.isfile(path):
        return params
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or ":" not in s:
                continue
            key, val = s.split(":", 1)
            val = re.sub(r"\s+#.*$", "", val).strip().strip("'\"")
            if val.lower() in ("true", "false"):
                val = val.lower() == "true"
            elif val.lower() in ("null", "~", ""):
                val = None
            params[key.strip()] = val
    return params


def read_run_script(path):
    """GENOME= / GTF= / SING_CACHE= / MASTER_SHEET= of the human branch of run_chipseq_slurm.sh."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    text = open(path).read()
    m = re.search(r'if \[ "\$\{SPECIES\}" = "Human" \]; then(.*?)(?:elif|else|fi)', text, re.S)
    block = m.group(1) if m else ""
    for key in ("GENOME", "GTF"):
        m = re.search(rf"^\s*{key}=['\"]([^'\"]+)['\"]", block, re.M)
        if m:
            out[key] = m.group(1)
    for key in ("SING_CACHE",):
        m = re.search(rf"^{key}=['\"]([^'\"]+)['\"]", text, re.M)
        if m:
            out[key] = m.group(1)
    return out


def scan_fasta(path, lengths=True):
    """[(name, description, length)] in file order. Uses a matching .fai when present (fast)."""
    fai = path + ".fai"
    if lengths and os.path.isfile(fai) and os.path.getmtime(fai) >= os.path.getmtime(path):
        seqs = []
        with open(fai) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                seqs.append((f[0], "", int(f[1])))
        return seqs, f"(lengths from {os.path.basename(fai)})"
    seqs, cur, bad = [], None, set()
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append(tuple(cur))
                h = line[1:].strip()
                parts = h.split(None, 1)
                cur = [parts[0] if parts else "", parts[1] if len(parts) > 1 else "", 0]
            elif cur is not None:
                s = line.strip()
                cur[2] += len(s)
                if lengths is not None and re.search(r"[^ACGTNacgtnRYKMSWBDHVrykmswbdhv*\-.]", s):
                    bad.add(cur[0])
    if cur:
        seqs.append(tuple(cur))
    note = f"(unexpected characters in: {short(bad)})" if bad else ""
    return seqs, note


def read_alias(path, column):
    """Same rules as the pipeline's rename_chroms.py: every name on a line maps to the name in 'column'."""
    header, rows = None, []
    with open_any(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").strip().split("\t") if "\t" in line else line.lstrip("#").split()
                continue
            rows.append([f.strip() for f in (line.rstrip("\n").split("\t") if "\t" in line else line.split())])
    if str(column).isdigit():
        idx = int(column) - 1
    elif header and column in header:
        idx = header.index(column)
    else:
        raise ValueError(f"column '{column}' not found in the alias table header {header}")
    mapping = {}
    for fields in rows:
        if idx >= len(fields) or not fields[idx]:
            continue
        for name in fields:
            if name:
                mapping[name] = fields[idx]
    return mapping


def apply_alias(names, mapping):
    renamed, unmapped, targets = [], [], {}
    for n in names:
        t = mapping.get(n, n)
        if n not in mapping and n not in set(mapping.values()):
            unmapped.append(n)
        renamed.append(t)
        targets.setdefault(t, []).append(n)
    collisions = {t: src for t, src in targets.items() if len(set(src)) > 1}
    return renamed, unmapped, collisions


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default="nextflow-params.yml", help="params file (default nextflow-params.yml)")
    ap.add_argument("--run-script", default="run_chipseq_slurm.sh", help="launch script to read GENOME/GTF from")
    ap.add_argument("--pipeline", default="chipseq", help="local pipeline checkout (default ./chipseq)")
    ap.add_argument("--fasta", help="host genome FASTA (default: human GENOME= in the run script)")
    ap.add_argument("--gtf", help="GTF (default: human GTF= in the run script)")
    ap.add_argument("--viral-fasta", help="viral FASTA (default: viral_fasta in the params file)")
    ap.add_argument("--refseq-fasta", help="RefSeq FASTA the GTF belongs to, to confirm identical chromosome lengths "
                                           "(default: the GCF_* FASTA next to the GTF, if present)")
    ap.add_argument("--skip-refseq-compare", action="store_true", help="do not read the RefSeq FASTA")
    a = ap.parse_args()

    params = read_params(a.params)
    run = read_run_script(a.run_script)
    fasta = a.fasta or run.get("GENOME") or params.get("fasta")
    gtf = a.gtf or run.get("GTF") or params.get("gtf")
    viral = a.viral_fasta or params.get("viral_fasta")
    alias_path = params.get("chrom_alias")
    alias_col = params.get("chrom_alias_column") or "ucsc"
    exclude = [x.strip() for x in str(params.get("exclude_contigs") or "").split(",") if x.strip()]
    blacklist = params.get("blacklist")
    bt2_index = params.get("bowtie2_index")

    print(f"Reference pre-flight on {os.uname().nodename} in {os.getcwd()}")
    print(f"  params file : {a.params if params else '(not found)'}")
    print(f"  host FASTA  : {fasta}\n  GTF         : {gtf}\n  viral FASTA : {viral}")
    print(f"  chrom_alias : {alias_path} (column '{alias_col}')\n  exclude     : {', '.join(exclude) or '-'}")

    # ---- 1. pipeline checkout and params ------------------------------------------------------------------------
    section("1. Pipeline checkout and params file")
    cfg = os.path.join(a.pipeline, "nextflow.config")
    known = set()
    if os.path.isfile(cfg):
        text = open(cfg).read()
        m = re.search(r"^params\s*\{(.*?)^\}", text, re.S | re.M)
        known = set(re.findall(r"^\s+(\w+)\s*=", m.group(1) if m else "", re.M))
        schema = os.path.join(a.pipeline, "nextflow_schema.json")
        if os.path.isfile(schema):      # some parameters (e.g. blacklist) are only declared in the schema
            import json

            def props(node):
                if isinstance(node, dict):
                    for k, v in node.items():
                        if k == "properties" and isinstance(v, dict):
                            known.update(v.keys())
                        props(v)
                elif isinstance(node, list):
                    for v in node:
                        props(v)
            props(json.load(open(schema)))
        need = ["modules/local/combine_viral_genome.nf", "modules/local/viral_read_stats.nf",
                "bin/combine_viral_genome.py", "modules/local/bowtie2_index_check.nf"]
        missing = [f for f in need if not os.path.isfile(os.path.join(a.pipeline, f))]
        if missing or "viral_fasta" not in known:
            say("FAIL", f"{a.pipeline}/ is an older pipeline version without --viral_fasta (missing: {short(missing) or 'viral_fasta param'}); "
                        "unpack the updated chipseq-2.1.0-modified.zip over it")
        else:
            say("PASS", f"{a.pipeline}/ has the viral-genome version of the pipeline")
        unknown = [k for k in params if known and k not in known]
        if unknown:
            say("WARN", f"params-file keys the pipeline does not know (typo?): {short(unknown)}")
        elif params:
            say("PASS", f"all {len(params)} params-file keys are pipeline parameters")
    else:
        say("FAIL", f"pipeline checkout not found: {cfg} (set --pipeline)")

    if not fasta or not gtf:
        say("FAIL", "host FASTA / GTF not known: give --fasta/--gtf or keep GENOME=/GTF= in the run script")
        return finish()

    # ---- 2. host FASTA -------------------------------------------------------------------------------------------
    section("2. Host genome FASTA")
    if not os.access(fasta, os.R_OK):
        say("FAIL", f"cannot read {fasta}")
        return finish()
    if os.path.basename(fasta).startswith("GCF_"):
        say("WARN", "this is the RefSeq (GCF_) FASTA, which has no mitochondrial sequence; use the GenBank (GCA_) "
                    "or UCSC hs1 FASTA to get chrM - and update GENOME= in run_chipseq_slurm.sh")
    print("      reading sequence names and lengths (a few minutes for a human genome) ...", flush=True)
    host, note = scan_fasta(fasta)
    names = [h[0] for h in host]
    lengths = {h[0]: h[2] for h in host}
    dups = sorted({n for n in names if names.count(n) > 1}) if len(names) < 5000 else []
    if dups:
        say("FAIL", f"duplicated host sequence names: {short(dups)}")
    say("PASS" if host else "FAIL", f"{len(host)} sequences, {sum(lengths.values()):,} bp {note}")
    mito = [h for h in host if h[2] == MITO_LEN or "mitochondr" in h[1].lower() or h[0] in ("chrM", "MT", "CP068254.1", "NC_012920.1")]
    if mito:
        say("PASS", f"mitochondrion present: {mito[0][0]} ({mito[0][2]:,} bp)")
    else:
        say("WARN" if not any(e in ("chrM", "MT") for e in exclude) else "FAIL",
            "no mitochondrial sequence (16,569 bp) in the host FASTA: mitochondrial reads cannot align to chrM")

    # ---- 3. chromosome alias ----------------------------------------------------------------------------------------
    section("3. Chromosome renaming (--chrom_alias)")
    mapping = {}
    final_host = names
    if alias_path:
        if not os.access(alias_path, os.R_OK):
            say("FAIL", f"cannot read alias table {alias_path}")
        else:
            try:
                mapping = read_alias(alias_path, alias_col)
            except ValueError as e:
                say("FAIL", str(e))
            if mapping:
                final_host, unmapped, coll = apply_alias(names, mapping)
                n_ren = sum(1 for n, f in zip(names, final_host) if n != f)
                if coll:
                    say("FAIL", f"different sequences would get the same name: {coll}")
                if n_ren == 0 and not set(names) & set(mapping.values()):
                    say("FAIL", "none of the host sequence names are in the alias table: wrong table or column?")
                else:
                    say("PASS", f"{n_ren} of {len(names)} sequences renamed, e.g. "
                                f"{', '.join(f'{n}->{f}' for n, f in list(zip(names, final_host))[:3])}")
                if unmapped:
                    say("WARN", f"{len(unmapped)} sequence(s) not in the alias table keep their names: {short(unmapped)}")
                if mito:
                    m_final = mapping.get(mito[0][0], mito[0][0])
                    say("PASS" if m_final == "chrM" else "WARN", f"mitochondrion will be named '{m_final}'")
    elif any(re.match(r"(CP|NC_|NT_|NW_|CM)\d", n) for n in names):
        say("WARN", "no --chrom_alias: sequences keep their accession names (e.g. CP068277.2) in every output")
    else:
        say("PASS", f"no --chrom_alias: sequence names used as they are ({short(names, 4)})")
    final_lengths = {f: lengths[n] for n, f in zip(names, final_host)}

    # ---- 4. same assembly as the RefSeq FASTA ------------------------------------------------------------------------
    section("4. Same assembly as the GTF's RefSeq genome")
    refseq = a.refseq_fasta
    if not refseq and gtf:
        cand = re.sub(r"\.gtf(\.gz)?$", ".fna", gtf)
        refseq = cand if os.path.isfile(cand) and os.path.abspath(cand) != os.path.abspath(fasta) else None
    if a.skip_refseq_compare or not refseq:
        say("WARN" if not a.skip_refseq_compare else "PASS",
            "RefSeq FASTA not compared (give --refseq-fasta to confirm identical chromosome lengths)"
            if not a.skip_refseq_compare else "comparison skipped (--skip-refseq-compare)")
    elif os.access(refseq, os.R_OK):
        print(f"      reading {os.path.basename(refseq)} ...", flush=True)
        ref, _ = scan_fasta(refseq)
        ref_final = {mapping.get(n, n): L for n, _, L in ref}
        diff = [c for c in ref_final if c in final_lengths and ref_final[c] != final_lengths[c]]
        absent = [c for c in ref_final if c not in final_lengths]
        if diff:
            say("FAIL", f"chromosome lengths differ from the RefSeq FASTA (different assembly?): {short(diff)}")
        elif absent:
            say("WARN", f"RefSeq sequences missing from the host FASTA after renaming: {short(absent)}")
        else:
            say("PASS", f"all {len(ref_final)} RefSeq chromosomes present with identical lengths "
                        f"(host FASTA adds: {short(set(final_lengths) - set(ref_final)) or 'nothing'})")
    else:
        say("WARN", f"cannot read {refseq}; lengths not compared")

    # ---- 5. GTF ------------------------------------------------------------------------------------------------------
    section("5. GTF annotation")
    if not os.access(gtf, os.R_OK):
        say("FAIL", f"cannot read {gtf}")
    else:
        maxend, nlines = {}, 0
        with open_any(gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 5)
                if len(f) < 5:
                    continue
                nlines += 1
                try:
                    e = int(f[4])
                except ValueError:
                    continue
                if e > maxend.get(f[0], 0):
                    maxend[f[0]] = e
        g_final = {mapping.get(c, c): e for c, e in maxend.items()}
        missing = [c for c in g_final if c not in final_lengths]
        outside = [c for c, e in g_final.items() if c in final_lengths and e > final_lengths[c]]
        matched = len(g_final) - len(missing)
        if matched == 0:
            hint = (" - with a pre-built index and no --chrom_alias, use the renamed GTF from genome/renamed/ as well"
                    if bt2_index and not alias_path else "")
            say("FAIL", f"no GTF chromosome matches the renamed genome (GTF: {short(maxend)}){hint}")
        else:
            say("PASS", f"{nlines:,} GTF lines on {len(g_final)} chromosomes; {matched} match the renamed genome")
        if missing:
            say("WARN", f"GTF chromosomes not in the genome (their genes are ignored): {short(missing)}")
        if outside:
            say("FAIL", f"genes extend beyond the chromosome end (GTF from another assembly?): {short(outside)}")
        elif matched:
            say("PASS", "all genes lie within the chromosome lengths")
        if mito and mapping.get(mito[0][0], mito[0][0]) not in g_final:
            print("INFO  the GTF has no mitochondrial genes (expected for the RefSeq T2T GTF): chrM peaks will be annotated as intergenic")

    # ---- 6. viral FASTA --------------------------------------------------------------------------------------------
    section("6. Viral genomes (--viral_fasta)")
    vnames = []
    if not viral:
        say("WARN", "no viral_fasta set: the run will align to the host genome only")
    elif not os.access(viral, os.R_OK):
        say("FAIL", f"cannot read {viral}")
    elif not re.search(r"\.(fa|fasta|fna)(\.gz)?$", viral):
        say("FAIL", f"{viral}: the pipeline accepts .fa/.fasta/.fna (optionally .gz) file names only")
    else:
        vseqs, vnote = scan_fasta(viral, lengths=True)
        vnames = [v[0] for v in vseqs]
        if not vseqs:
            say("FAIL", f"{viral} contains no FASTA sequences")
        else:
            for n, d, L in vseqs:
                print(f"      {n:<20} {L:>9,} bp  {d[:70]}")
            dup = sorted({n for n in vnames if vnames.count(n) > 1})
            empty = [n for n, _, L in vseqs if L == 0]
            clash = [n for n in vnames if n in set(final_host) or n in set(names)]
            odd = [n for n in vnames if not re.fullmatch(r"[A-Za-z0-9_.|+\-]+", n)]
            small = [n for n, _, L in vseqs if 0 < L < 1000]
            if dup:
                say("FAIL", f"duplicated viral sequence names (first word of the header): {short(dup)}")
            if empty:
                say("FAIL", f"empty viral sequences: {short(empty)}")
            if clash:
                say("FAIL", f"viral names already used by host sequences: {short(clash)}")
            if vnote:
                say("FAIL", f"viral FASTA {vnote}")
            if not (dup or empty or clash or vnote):
                say("PASS", f"{len(vseqs)} viral sequence(s), {sum(v[2] for v in vseqs):,} bp, unique names, no clash with host names")
            if odd:
                say("WARN", f"unusual characters in viral names (prefer letters, digits, . _ -): {short(odd)}")
            if small:
                say("WARN", f"very short viral sequences (< 1 kb; partial genome or wrong file?): {short(small)}")
            rna = [n for n, d, _ in vseqs if re.search(r"influenza|zika|dengue|coronavirus|sars|measles|vesicular", d, re.I)]
            if rna:
                print(f"INFO  {short(rna)} look like RNA-virus genomes: DNA ChIP should contain (almost) no reads from them")
            total = sum(final_lengths.values()) + sum(v[2] for v in vseqs)
            say("PASS" if total < 4_294_967_096 else "WARN",
                f"combined genome {total:,} bp: {'standard' if total < 4_294_967_096 else 'LARGE (.bt2l)'} Bowtie2 index")

    # ---- 7. exclude_contigs ---------------------------------------------------------------------------------------------
    section("7. Sequences to exclude (--exclude_contigs)")
    if not exclude:
        print("INFO  none set")
    else:
        allnames = set(final_host) | set(vnames)
        for e in exclude:
            say("PASS" if e in allnames else "WARN",
                f"'{e}' found" if e in allnames else
                f"'{e}' is not a sequence of the final genome, so the pipeline will ignore it "
                f"(names are the final, renamed ones: {short(sorted(allnames), 6)})")

    # ---- 8. blacklist ------------------------------------------------------------------------------------------------------
    section("8. Blacklist")
    if not blacklist:
        print("INFO  no blacklist set")
    elif not os.access(blacklist, os.R_OK):
        say("FAIL", f"cannot read {blacklist}")
    else:
        bad, chroms, n = [], set(), 0
        with open_any(blacklist) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) < 3 or line.startswith(("#", "track", "browser")):
                    continue
                n += 1
                c = mapping.get(f[0], f[0])
                chroms.add(c)
                if c in final_lengths and int(f[2]) > final_lengths[c]:
                    bad.append(c)
        miss = [c for c in chroms if c not in final_lengths]
        say("FAIL" if chroms and not set(chroms) & set(final_lengths) else "PASS",
            f"{n} intervals on {len(chroms)} chromosomes; {len(chroms) - len(miss)} match the genome")
        if miss:
            say("WARN", f"blacklist chromosomes not in the genome (ignored): {short(miss)}")
        if bad:
            say("FAIL", f"blacklist intervals beyond chromosome ends (other assembly?): {short(set(bad))}")

    # ---- 9. pre-built Bowtie2 index -----------------------------------------------------------------------------------
    section("9. Pre-built Bowtie2 index")
    if not bt2_index:
        print("INFO  no --bowtie2_index: the combined index is built during the run (about 1-2 h for human)")
    else:
        if alias_path:
            say("FAIL", "--bowtie2_index cannot be combined with --chrom_alias (the pipeline stops): pass the renamed "
                        "FASTA from genome/renamed/ as --fasta and drop chrom_alias")
        inspect = shutil.which("bowtie2-inspect")
        if not inspect:
            say("WARN", "bowtie2-inspect not on PATH: index contents not checked here (the pipeline checks them before aligning)")
        else:
            idx = None
            for root, _, files in os.walk(bt2_index):
                for f in files:
                    if f.endswith((".rev.1.bt2", ".rev.1.bt2l")):
                        idx = os.path.join(root, re.sub(r"\.rev\.1\.bt2l?$", "", f))
            if not idx:
                say("FAIL", f"no Bowtie2 index files in {bt2_index}")
            else:
                out = subprocess.run([inspect, "-n", idx], capture_output=True, text=True)
                inames = {l.split()[0] for l in out.stdout.splitlines() if l.strip()}
                need = set(final_host) | set(vnames)
                lack = need - inames
                say("FAIL" if lack else "PASS",
                    f"index lacks {len(lack)} sequence(s): {short(sorted(lack))}" if lack else
                    f"index contains all {len(need)} host + viral sequences")

    # ---- 10. environment -------------------------------------------------------------------------------------------------
    section("10. Environment")
    pre = "prefetch_singularity_images.sh"
    if os.path.isfile(pre):
        env = dict(os.environ, PIPELINE_DIR=a.pipeline)
        if run.get("SING_CACHE"):
            env["CACHE_DIR"] = run["SING_CACHE"]
        r = subprocess.run(["bash", pre, "--check"], capture_output=True, text=True, env=env)
        say("PASS" if r.returncode == 0 else "WARN",
            r.stdout.strip().splitlines()[0] if r.stdout.strip() else (r.stderr.strip() or "image check failed"))
    else:
        print(f"INFO  {pre} not found: Singularity image cache not checked")
    free = shutil.disk_usage(".").free / 1e9
    say("PASS" if free > 200 else "WARN", f"{free:,.0f} GB free in the launch directory "
        f"(a human combined index needs ~4 GB; each batch's work directory can reach 100+ GB)")
    return finish()


def finish():
    print(f"\nSUMMARY: {RESULTS['PASS']} passed, {RESULTS['WARN']} warnings, {RESULTS['FAIL']} failed")
    print("=> Ready to run." if RESULTS["FAIL"] == 0 else "=> Fix the FAIL lines before running the pipeline.")
    sys.exit(1 if RESULTS["FAIL"] else 0)


if __name__ == "__main__":
    main()
