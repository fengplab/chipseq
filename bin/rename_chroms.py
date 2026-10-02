#!/usr/bin/env python3

"""
Rename chromosome/sequence names in a FASTA, GTF/GFF or BED file using a chromosome alias table,
e.g. RefSeq 'NC_000001.11' / 'NC_060925.1' -> UCSC 'chr1'.

Alias table: every name on a line is an equivalent name for the same sequence (UCSC <assembly>.chromAlias.txt
format; tab- or whitespace-separated). The name to convert TO is taken from --column: either a column name
from a '#'-prefixed header line (e.g. 'ucsc') or a 1-based column number.

Names that are not in the alias table are kept unchanged and reported. The script fails if no sequence was
renamed (usually the wrong alias table or column) or if two different input names would map to the same
output name.
"""

import argparse
import gzip
import sys


def open_any(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "r")


def read_alias(path, column):
    header = None
    rows = []
    with open_any(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").strip().split("\t") if "\t" in line else line.lstrip("#").split()
                continue
            fields = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            rows.append([f.strip() for f in fields])

    if column.isdigit():
        idx = int(column) - 1
    elif header is not None and column in [h.strip() for h in header]:
        idx = [h.strip() for h in header].index(column)
    elif header is None and column == "ucsc":
        idx = 0
    else:
        sys.exit(
            "ERROR: column '{}' not found in the alias table header {}. Use --chrom_alias_column with a column "
            "name or a 1-based column number.".format(column, header)
        )

    mapping = {}
    for fields in rows:
        if idx >= len(fields) or not fields[idx]:
            continue
        target = fields[idx]
        for name in fields:
            if name:
                mapping[name] = target
    return mapping, header, idx


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="FASTA, GTF/GFF or BED file (optionally gzipped)")
    parser.add_argument("alias", help="Chromosome alias table")
    parser.add_argument("output", help="Output file (uncompressed)")
    parser.add_argument("--format", choices=["fasta", "tab"], required=True)
    parser.add_argument("--column", default="ucsc", help="Alias column to rename to (name or 1-based number)")
    parser.add_argument("--report", default=None, help="Write a per-name renaming report here")
    parser.add_argument("--allow_none", action="store_true", help="Do not fail if nothing was renamed")
    a = parser.parse_args(args)

    mapping, header, idx = read_alias(a.alias, a.column)
    sys.stderr.write("Alias table: {} names -> column {} ({})\n".format(
        len(mapping), idx + 1, header[idx] if header and idx < len(header) else "no header"))

    seen = {}          # original name -> new name
    reverse = {}       # new name -> original name (to detect collisions)
    counts = {}

    def rename(name):
        new = mapping.get(name, name)
        if name not in seen:
            seen[name] = new
            other = reverse.get(new)
            if other is not None and other != name:
                sys.exit("ERROR: '{}' and '{}' would both be renamed to '{}'.".format(other, name, new))
            reverse[new] = name
        counts[name] = counts.get(name, 0) + 1
        return new

    with open_any(a.input) as fin, open(a.output, "w") as fout:
        if a.format == "fasta":
            for line in fin:
                if line.startswith(">"):
                    parts = line[1:].rstrip("\n").split(None, 1)
                    name = parts[0] if parts else ""
                    rest = (" " + parts[1]) if len(parts) > 1 else ""
                    fout.write(">" + rename(name) + rest + "\n")
                else:
                    fout.write(line)
        else:
            for line in fin:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    fout.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                fields[0] = rename(fields[0])
                fout.write("\t".join(fields) + "\n")

    renamed = sorted(n for n, t in seen.items() if t != n)
    kept = sorted(n for n, t in seen.items() if t == n)
    if a.report:
        with open(a.report, "w") as rep:
            rep.write("original_name\tnew_name\tstatus\tn_records\n")
            for n in renamed:
                rep.write("{}\t{}\trenamed\t{}\n".format(n, seen[n], counts[n]))
            for n in kept:
                status = "already_target" if n in mapping else "not_in_alias_table"
                rep.write("{}\t{}\t{}\t{}\n".format(n, n, status, counts[n]))

    unaliased = [n for n in kept if n not in mapping]
    sys.stderr.write("{} names renamed, {} already in target form, {} not in alias table (kept)\n".format(
        len(renamed), len(kept) - len(unaliased), len(unaliased)))
    if unaliased:
        sys.stderr.write("  not in alias table, e.g.: {}\n".format(" ".join(unaliased[:5])))
    if seen and not renamed and not (set(seen) & set(mapping)):
        if a.allow_none:
            sys.stderr.write(
                "WARNING: none of the {} sequence names ({}) are in the alias table; original names were kept. "
                "Check that the alias table belongs to this assembly.\n".format(len(seen), " ".join(sorted(seen)[:5])))
            return 0
        sys.exit(
            "ERROR: none of the {} sequence names ({}) are in the alias table. Check that the alias table belongs "
            "to this assembly.".format(len(seen), " ".join(sorted(seen)[:5]))
        )


if __name__ == "__main__":
    sys.exit(main())
