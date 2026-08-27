#!/usr/bin/env python3
"""
Build gene x sample count and RPKM matrices from Flow `*.summary_gene.tsv` files.

Each input has columns: "Gene name (Gene ID)", Length, cDNA #, cDNA %. Gene length is
supplied per row, so RPKM needs no GTF:

    RPKM = cDNA# * 1e9 / (length_bp * library_total)

`library_total` is the sum of EVERY row's cDNA#, intergenic included, because that is the
library the percentages are computed against - verified: summing the column reproduces the
total implied by `cDNA %` exactly. Intergenic is counted in that denominator but is not a
gene, so it is not emitted as a matrix row.

Rows whose symbol is "None" (unnamed genes) are dropped. Symbols mapping to several Ensembl
IDs are collapsed by summing counts and taking the longest annotation; only 57 symbols are
affected and they carry ~4% of the library, most of it small-RNA families (Y_RNA, U3,
Metazoa_SRP) that no symbol-keyed gene list will ask for.

Samples differ in which genes they report at all, so the matrix is the union across files
with 0 filled in for absent genes.

Python 3.6 compatible.
"""

import argparse
import os
import re
import sys

PAT = re.compile(r"^(.*)\s+\(([^)]*)\)$")
ORDER_HINT = ["Centro", "NT", "No"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", "--dir", required=True, help="Directory of *.summary_gene.tsv files")
    p.add_argument("-o", "--outdir", required=True, help="Where to write counts/rpkm/lengths TSVs")
    p.add_argument("--prefix", default="centro", help="Output filename prefix (default centro)")
    return p.parse_args()


def sample_name(fn):
    return os.path.basename(fn).replace(".summary_gene.tsv", "")


def sort_key(name):
    for i, h in enumerate(ORDER_HINT):
        if name.startswith(h):
            return (i, name)
    return (len(ORDER_HINT), name)


def read_one(path):
    """Return (counts_by_symbol, lengths_by_symbol, library_total)."""
    counts, lengths, total = {}, {}, 0
    with open(path) as fh:
        header = fh.readline()
        if "cDNA" not in header:
            sys.exit("%s does not look like a summary_gene.tsv (header: %r)" % (path, header[:60]))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                c = int(f[2])
                L = int(f[1])
            except ValueError:
                continue
            total += c                      # library denominator includes intergenic
            m = PAT.match(f[0])
            if not m:
                continue
            sym = m.group(1).strip()
            if sym in ("intergenic", "None", ""):
                continue
            counts[sym] = counts.get(sym, 0) + c
            lengths[sym] = max(lengths.get(sym, 0), L)
    return counts, lengths, total


def main():
    args = parse_args()
    files = sorted([os.path.join(args.dir, f) for f in os.listdir(args.dir)
                    if f.endswith(".summary_gene.tsv")], key=lambda p: sort_key(sample_name(p)))
    if not files:
        sys.exit("no *.summary_gene.tsv under %s" % args.dir)

    per, lengths, totals = {}, {}, {}
    for path in files:
        s = sample_name(path)
        c, L, t = read_one(path)
        per[s], totals[s] = c, t
        for g, v in L.items():
            lengths[g] = max(lengths.get(g, 0), v)
        print("  %-22s genes=%6d  library total=%12s" % (s, len(c), "{:,}".format(t)))

    samples = [sample_name(p) for p in files]
    genes = sorted(set(g for c in per.values() for g in c))
    print("\nunion across samples: %d genes x %d samples" % (len(genes), len(samples)))

    os.makedirs(args.outdir, exist_ok=True)
    cpath = os.path.join(args.outdir, args.prefix + "_gene_counts.tsv")
    rpath = os.path.join(args.outdir, args.prefix + "_gene_rpkm.tsv")
    with open(cpath, "w") as fc, open(rpath, "w") as fr:
        fc.write("gene\t" + "\t".join(samples) + "\n")
        fr.write("gene\t" + "\t".join(samples) + "\n")
        for g in genes:
            L = lengths.get(g, 0)
            cs = [per[s].get(g, 0) for s in samples]
            fc.write(g + "\t" + "\t".join(str(v) for v in cs) + "\n")
            if L > 0:
                rk = [v * 1e9 / (L * float(totals[s])) if totals[s] else 0.0
                      for v, s in zip(cs, samples)]
                fr.write(g + "\t" + "\t".join("%.6g" % v for v in rk) + "\n")
    nolen = sum(1 for g in genes if lengths.get(g, 0) <= 0)
    print("wrote %s" % cpath)
    print("wrote %s%s" % (rpath, "" if not nolen else "  (%d genes had no length and were omitted)" % nolen))


if __name__ == "__main__":
    main()
