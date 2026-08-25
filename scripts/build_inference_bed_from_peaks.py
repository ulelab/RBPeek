#!/usr/bin/env python3
"""
Build an inference BED (`binf`) for intersect_inference_bed.py from replicate Clippy peak
calls. Generalises the THRAP3-specific build_thrap3_inference_bed.py to any bait.

Three transformations are applied, each for a reason that would otherwise break the run:

1. Chromosome naming. Clippy run against an Ensembl FASTA emits "1"; the RBPeek panel and
   decoys are UCSC-style "chr1". intersect_inference_bed.py normalises names only on its
   internal merge path, never for the inference BED, so an unnormalised BED intersects
   nothing and the run completes silently with an empty result.

2. Replicate reproducibility. Replicates usually differ several-fold in depth, so a plain
   union is dominated by whichever library was sequenced deepest. Peaks are merged
   strand-aware and kept only where at least --min-reps replicates contribute.

3. Single-nucleotide anchors. intersect_inference_bed.py assigns each panel interval's whole
   score to one offset relative to the locus, so an inference locus wider than 1 bp smears
   the metaprofile by its own width. Each merged region collapses to its midpoint.

Output BED6: chrom, start, end(=start+1), <label>_<n>reps_<i>, <n>, strand
The score column carries replicate support so results can be split by reproducibility tier.

Written to run under Python 3.6 (the HPC login node's system python), so no f-string debug
specifiers, no walrus, no builtin generics, no subprocess text=/capture_output=.
"""

import argparse
import subprocess
import sys
from collections import Counter
from pathlib import Path

PRIMARY = set(["chr" + str(c) for c in range(1, 23)] + ["chrX", "chrY", "chrM"])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", "--peaks", nargs="+", required=True, type=Path,
                   help="Replicate peak BEDs, in replicate order")
    p.add_argument("-l", "--label", required=True,
                   help="Bait label used in output filenames and locus names (e.g. Centro_apex)")
    p.add_argument("--min-reps", type=int, default=2,
                   help="Keep merged regions supported by at least this many replicates (default 2)")
    p.add_argument("--anchor", choices=["midpoint", "start"], default="midpoint",
                   help="How to collapse each merged region to 1 nt (default midpoint)")
    p.add_argument("--keep-scaffolds", action="store_true",
                   help="Keep non-primary contigs; by default only chr1-22,X,Y,M are kept")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory")
    return p.parse_args()


def require(tool):
    from shutil import which
    if which(tool) is None:
        sys.exit("%s not found on PATH. Try: conda activate rbpeek" % tool)


def normalise(peaks, workdir, keep_scaffolds):
    """Rewrite each replicate as sorted, chr-prefixed BED6 tagged with its replicate label."""
    out = []
    for idx, src in enumerate(peaks, 1):
        if not src.is_file():
            sys.exit("peak file not found: %s" % src)
        rep = "R%d" % idx
        dst = workdir / ("%s.bed" % rep)
        unsorted = workdir / ("%s.unsorted" % rep)
        kept = dropped = 0
        with open(str(src)) as fin, open(str(unsorted), "w") as fout:
            for line in fin:
                if not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 6:
                    sys.exit("%s: expected >=6 columns, got %d" % (src.name, len(c)))
                chrom = c[0] if c[0].startswith("chr") else "chr" + c[0]
                if chrom == "chrMT":
                    chrom = "chrM"
                if not keep_scaffolds and chrom not in PRIMARY:
                    dropped += 1
                    continue
                kept += 1
                # col4 carries the replicate tag; bedtools merge -o distinct on it is what
                # produces the per-region replicate support count.
                fout.write("\t".join([chrom, c[1], c[2], rep, c[4], c[5]]) + "\n")
        with open(str(dst), "w") as fout:
            subprocess.check_call(["sort", "-k1,1", "-k2,2n", str(unsorted)], stdout=fout)
        unsorted.unlink()
        print("  %-4s %-52s kept=%6d  dropped_scaffold=%d" % (rep, src.name[:52], kept, dropped))
        out.append(dst)
    return out


def main():
    args = parse_args()
    require("bedtools")
    require("sort")
    if args.min_reps < 1 or args.min_reps > len(args.peaks):
        sys.exit("--min-reps must be between 1 and the number of peak files (%d)" % len(args.peaks))

    outdir = args.outdir
    workdir = outdir / "work"
    workdir.mkdir(parents=True, exist_ok=True)
    n_reps = len(args.peaks)

    print("[1/4] normalising chromosome names (%d replicates)" % n_reps)
    reps = normalise(args.peaks, workdir, args.keep_scaffolds)

    print("[2/4] merging strand-aware across replicates")
    allbed = workdir / "all.bed"
    with open(str(allbed), "w") as fout:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n"] + [str(p) for p in reps], stdout=fout)
    merged = workdir / "merged_raw.bed"
    with open(str(merged), "w") as fout:
        subprocess.check_call(["bedtools", "merge", "-i", str(allbed), "-s",
                               "-c", "4,5,6", "-o", "distinct,sum,distinct"], stdout=fout)

    print("[3/4] applying reproducibility filter and collapsing to 1 nt anchors")
    hist = Counter()
    combos = Counter()
    widths = []
    rows = []
    with open(str(merged)) as fin:
        for line in fin:
            c = line.rstrip("\n").split("\t")
            chrom, start, end = c[0], int(c[1]), int(c[2])
            reps_here = sorted(set(c[3].split(",")))
            strand = c[5].split(",")[0]
            n = len(reps_here)
            hist[n] += 1
            if n < args.min_reps:
                continue
            combos[",".join(reps_here)] += 1
            widths.append(end - start)
            anchor = start if args.anchor == "start" else (start + end) // 2
            rows.append((chrom, anchor, strand, n))

    rows.sort(key=lambda r: (r[0], r[1]))
    # Midpoint collapse can land two nearby regions on the same anchor; keep loci unique.
    seen = set()
    final = []
    for chrom, anchor, strand, n in rows:
        key = (chrom, anchor, strand)
        if key in seen:
            continue
        seen.add(key)
        final.append((chrom, anchor, strand, n))

    out_bed = outdir / ("%s_merged_min%drep_anchors.bed" % (args.label, args.min_reps))
    with open(str(out_bed), "w") as fout:
        for i, (chrom, anchor, strand, n) in enumerate(final, 1):
            fout.write("\t".join([chrom, str(anchor), str(anchor + 1),
                                  "%s_%dreps_%d" % (args.label, n, i), str(n), strand]) + "\n")

    print("[4/4] summary")
    total = sum(hist.values())
    print("  merged regions (>=1 rep): %d" % total)
    for k in sorted(hist):
        print("    present in %d/%d reps: %7d (%5.1f%%)" % (k, n_reps, hist[k], 100.0 * hist[k] / total))
    cum = sum(v for k, v in hist.items() if k >= args.min_reps)
    print("  passing >=%d reps: %d (%.1f%%)" % (args.min_reps, cum, 100.0 * cum / total))
    print("  replicate combinations in the kept set:")
    for combo, n in combos.most_common():
        print("    %-16s %6d" % (combo, n))
    if widths:
        print("  merged width: mean=%.1f max=%d" % (sum(widths) / float(len(widths)), max(widths)))
    print("  anchors written: %d (deduped from %d)" % (len(final), len(rows)))
    print("  -> %s" % out_bed)


if __name__ == "__main__":
    main()
