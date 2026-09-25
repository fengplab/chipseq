#!/usr/bin/env python3

"""
Convert RepeatMasker / CenSat annotations to a sorted, chromosome-size-aware BED6+2 file
(chrom, start, end, name, score, strand, class, family) ready for bedToBigBed with
assets/feature_bed6plus2.as.

Supported inputs (optionally gzipped, auto-detected with --format auto):
  rmsk_out   RepeatMasker *.out file (3 header lines, whitespace separated)
  ucsc_rmsk  UCSC rmsk table dump (rmsk.txt[.gz]; 17 columns starting with 'bin')
  bed        Any BED3+ file (e.g. the T2T CenSat BED). Class/family are taken from columns 7/8 when
             they look like text (BED6+2 produced by this script), otherwise derived from the name.
"""

import argparse
import gzip
import re
import sys


def open_any(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1].isdigit():
                sizes[fields[0]] = int(fields[1])
    return sizes


def censat_class(name):
    """'hor_1(S1C1/5/19H1L)' -> 'hor'; 'active_hor(S3C1H1L)' -> 'active_hor'; 'hsat2_2' -> 'hsat2'; 'ct_1_1(p_arm)' -> 'ct'."""
    base = re.sub(r"\(.*$", "", name)
    base = re.sub(r"(_\d+)+$", "", base)
    return base if base else name


def split_repeat_name(name):
    """RepeatMasker style 'AluY#SINE/Alu' -> ('AluY', 'SINE', 'Alu')."""
    if "#" in name:
        rep, cls = name.split("#", 1)
        fam = cls.split("/", 1)[1] if "/" in cls else cls
        return rep, cls.split("/", 1)[0], fam
    return name, None, None


def split_class_family(class_family):
    """'SINE/Alu' -> ('SINE', 'Alu'); 'Simple_repeat' -> ('Simple_repeat', 'Simple_repeat')."""
    if "/" in class_family:
        cls, fam = class_family.split("/", 1)
        return cls, fam
    return class_family, class_family


def detect_format(path):
    with open_any(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            s = line.strip()
            if s.startswith("SW") or s.startswith("score") or s.startswith("There were no repetitive"):
                return "rmsk_out"
            if s.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 17 and fields[0].isdigit() and fields[6].isdigit() and fields[7].isdigit():
                return "ucsc_rmsk"
            ws = s.split()
            if len(ws) >= 15 and ws[5].isdigit() and ws[6].isdigit() and ws[8] in ("+", "C"):
                return "rmsk_out"
            return "bed"
    return "bed"


def parse_rmsk_out(fh):
    for line in fh:
        ws = line.split()
        if len(ws) < 15 or not ws[0].isdigit():
            continue
        chrom, start, end = ws[4], int(ws[5]) - 1, int(ws[6])
        strand = "-" if ws[8] == "C" else "+"
        rep_name = ws[9]
        cls, fam = split_class_family(ws[10])
        score = int(ws[0])
        yield chrom, start, end, rep_name, score, strand, cls, fam


def parse_ucsc_rmsk(fh):
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 13 or not f[6].isdigit():
            continue
        yield f[5], int(f[6]), int(f[7]), f[10], int(f[1]), f[9], f[11], f[12]


def parse_bed(fh, feature_type):
    for line in fh:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3 or not (f[1].isdigit() and f[2].isdigit()):
            continue
        chrom, start, end = f[0], int(f[1]), int(f[2])
        name = f[3] if len(f) > 3 and f[3] else "{}:{}-{}".format(chrom, start, end)
        score = f[4] if len(f) > 4 else "0"
        strand = f[5] if len(f) > 5 and f[5] in ("+", "-", ".") else "."
        cls = fam = None
        # BED6+2 written by this script (or compatible): text columns 7 and 8
        if len(f) >= 8 and not re.match(r"^-?\d+$", f[6]) and not re.match(r"^-?\d+$", f[7]):
            cls, fam = f[6], f[7]
        if cls is None:
            rep, c, fm = split_repeat_name(name)
            if c is not None:
                name, cls, fam = rep, c, fm
        if cls is None:
            cls = censat_class(name) if feature_type == "censat" else name
            fam = name
        yield chrom, start, end, name, score, strand, cls, fam


def clean_score(score):
    try:
        return max(0, min(1000, int(float(score))))
    except (TypeError, ValueError):
        return 0


def clean_text(text):
    text = str(text).strip().replace("\t", "_").replace(" ", "_")
    return text if text else "NA"


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="RepeatMasker .out, UCSC rmsk.txt or BED file (optionally gzipped)")
    parser.add_argument("chrom_sizes", help="Chromosome sizes file (chrom<TAB>size)")
    parser.add_argument("output", help="Output BED6+2 file")
    parser.add_argument("--format", choices=["auto", "rmsk_out", "ucsc_rmsk", "bed"], default="auto")
    parser.add_argument("--type", choices=["repeatmasker", "censat", "generic"], default="generic")
    a = parser.parse_args(args)

    sizes = read_chrom_sizes(a.chrom_sizes)
    fmt = detect_format(a.input) if a.format == "auto" else a.format
    sys.stderr.write("Input format: {}\n".format(fmt))

    records, skipped_chrom, skipped_len = [], 0, 0
    with open_any(a.input) as fh:
        if fmt == "rmsk_out":
            it = parse_rmsk_out(fh)
        elif fmt == "ucsc_rmsk":
            it = parse_ucsc_rmsk(fh)
        else:
            it = parse_bed(fh, a.type)
        for chrom, start, end, name, score, strand, cls, fam in it:
            if chrom not in sizes:
                skipped_chrom += 1
                continue
            start, end = max(0, start), min(end, sizes[chrom])
            if end <= start:
                skipped_len += 1
                continue
            records.append((chrom, start, end, clean_text(name), clean_score(score), strand, clean_text(cls), clean_text(fam)))

    # bedToBigBed requires 'sort -k1,1 -k2,2n' ordering (byte-wise chromosome names)
    records.sort(key=lambda r: (r[0].encode(), r[1], r[2]))
    with open(a.output, "w") as out:
        for r in records:
            out.write("\t".join(map(str, r)) + "\n")

    sys.stderr.write(
        "Wrote {} features; skipped {} on chromosomes absent from the genome and {} empty after clipping\n".format(
            len(records), skipped_chrom, skipped_len
        )
    )
    if not records:
        sys.exit("ERROR: no features left after matching to the genome chromosome names. Check chromosome naming (e.g. 'chr1' vs '1').")


if __name__ == "__main__":
    sys.exit(main())
