#!/usr/bin/env python3

"""
Append viral genome(s) to a host genome FASTA so that reads can be aligned to both at once.

  - Host sequences are copied unchanged (after any --chrom_alias renaming done upstream).
  - Viral sequence names are cut to the first word of each FASTA header (the rest is kept as the description
    in the viral contig table), Windows line endings are removed and sequence lines are kept as given.
  - The run stops if the viral FASTA has no sequences, an empty sequence, a duplicated name, or a name that is
    already used by a host sequence (reads could then not be assigned to host or virus).

Outputs:
  <out>                       combined FASTA (host first, then viral sequences)
  viral_contigs.tsv           name, length, description of every viral sequence (no header line)
"""

import argparse
import gzip
import re
import sys


def open_any(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "r")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", required=True)
    ap.add_argument("--viral", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contigs", default="viral_contigs.tsv")
    a = ap.parse_args()

    # ---- viral sequences -------------------------------------------------------------
    viral = []          # [name, description, [lines], length]
    with open_any(a.viral) as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    sys.exit(f"ERROR: {a.viral}: a FASTA header has no sequence name.")
                parts = header.split(None, 1)
                viral.append([parts[0], parts[1] if len(parts) > 1 else "", [], 0])
            else:
                if not viral:
                    sys.exit(f"ERROR: {a.viral}: sequence data found before the first '>' header - not a FASTA file?")
                seq = line.strip()
                if re.search(r"[^A-Za-z*\-.]", seq):
                    sys.exit(f"ERROR: {a.viral}: unexpected characters in the sequence of '{viral[-1][0]}'.")
                viral[-1][2].append(seq)
                viral[-1][3] += len(seq)
    if not viral:
        sys.exit(f"ERROR: {a.viral} contains no FASTA sequences.")
    empty = [v[0] for v in viral if v[3] == 0]
    if empty:
        sys.exit(f"ERROR: {a.viral}: empty sequence(s): {', '.join(empty)}")
    seen, dup = set(), []
    for v in viral:
        (dup.append(v[0]) if v[0] in seen else seen.add(v[0]))
    if dup:
        sys.exit(f"ERROR: {a.viral}: duplicated sequence name(s): {', '.join(sorted(set(dup)))}. "
                 "Names are the first word of each header and must be unique.")
    odd = [v[0] for v in viral if not re.fullmatch(r"[A-Za-z0-9_.|+\-]+", v[0])]
    if odd:
        sys.stderr.write("WARNING: unusual characters in viral sequence name(s) "
                         f"{', '.join(odd)}; some tools may not accept them.\n")

    # ---- host: copy unchanged, collecting names ------------------------------------------------
    host_names = set()
    n_host = 0
    with open_any(a.host) as fin, open(a.out, "w") as fout:
        last = "\n"
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                host_names.add(name)
                n_host += 1
            fout.write(line)
            last = line
        if not last.endswith("\n"):
            fout.write("\n")
        clash = [v[0] for v in viral if v[0] in host_names]
        if clash:
            fout.close()
            sys.exit(f"ERROR: viral sequence name(s) {', '.join(clash)} already exist in the host genome {a.host}. "
                     "Rename them in the viral FASTA header.")
        for name, desc, lines, _ in viral:
            fout.write(f">{name}\n")
            for seq in lines:
                fout.write(seq + "\n")

    with open(a.contigs, "w") as fh:
        for name, desc, _, length in viral:
            fh.write(f"{name}\t{length}\t{desc}\n")

    total = sum(v[3] for v in viral)
    sys.stderr.write(f"Added {len(viral)} viral sequence(s) ({total:,} bp) to {n_host} host sequence(s): "
                     f"{', '.join(v[0] for v in viral)}\n")


if __name__ == "__main__":
    main()
