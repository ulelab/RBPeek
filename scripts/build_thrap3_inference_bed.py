#!/usr/bin/env python3
"""
Build the THRAP3 inference BED for intersect_inference_bed.py from replicate Clippy peak calls.

  1. Chromosome names are converted to UCSC style (chr-prefixed; MT -> chrM) and, unless
     --keep-scaffolds is given, records on non-primary contigs are removed.
  2. Peaks from all replicates are pooled and overlapping peaks on the same strand are merged.
  3. Merged regions supported by fewer than --min-reps distinct replicates are discarded.
  4. Each remaining region is reduced to a 1 nt anchor, by default its midpoint.

Output BED6:
    chrom  start  end (= start + 1)  name  score  strand
where name is THRAP3_<n>reps_<i> and score is the number of supporting replicates.
"""

import argparse
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_RAW = REPO / "THRAP3" / "raw"
# Sample name -> replicate label used in the peak file names.
SAMPLE_TO_REP = {
    "THRAP3_1": "R1",
    "THRAP3_2": "R2",
    "THRAP3_3": "R3",
    "THRAP3_L": "R4",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW,
                   help="Directory containing the replicate *_genome.*_Peaks.bed files")
    p.add_argument("--min-reps", type=int, default=2,
                   help="Keep merged regions supported by at least this many replicates (default 2)")
    p.add_argument("--anchor", choices=["midpoint", "start"], default="midpoint",
                   help="How to collapse each merged region to 1 nt (default midpoint)")
    p.add_argument("--keep-scaffolds", action="store_true",
                   help="Keep non-primary contigs; by default only chr1-22,X,Y,M are kept")
    p.add_argument("-o", "--outdir", type=Path, default=REPO / "THRAP3",
                   help="Output directory (default THRAP3/)")
    return p.parse_args()


PRIMARY = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M", "MT"]}


def normalise(raw_dir, workdir, keep_scaffolds):
    """Rewrite each replicate's peaks as sorted, chr-prefixed BED6 tagged with its replicate."""
    out = []
    for sample, rep in sorted(SAMPLE_TO_REP.items(), key=lambda kv: kv[1]):
        hits = sorted(raw_dir.glob(f"THRAP3_{rep}_genome.*_Peaks.bed"))
        if len(hits) != 1:
            sys.exit(f"expected exactly one genome Peaks.bed for {rep} in {raw_dir}, found {len(hits)}")
        src = hits[0]
        dst = workdir / f"{rep}.chr.bed"
        kept = dropped = 0
        unsorted = workdir / f"{rep}.unsorted"
        with open(src) as fin, open(unsorted, "w") as fout:
            for line in fin:
                if not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 6:
                    sys.exit(f"{src.name}: expected >=6 columns, got {len(c)}")
                chrom = c[0] if c[0].startswith("chr") else "chr" + c[0]
                if chrom == "chrMT":
                    chrom = "chrM"
                if not keep_scaffolds and chrom not in PRIMARY:
                    dropped += 1
                    continue
                kept += 1
                # Column 4 holds the replicate label, so that bedtools merge -o distinct
                # reports the replicates supporting each merged region.
                fout.write("\t".join([chrom, c[1], c[2], rep, c[4], c[5]]) + "\n")
        with open(dst, "w") as fout:
            subprocess.run(["sort", "-k1,1", "-k2,2n", str(unsorted)], stdout=fout, check=True)
        unsorted.unlink()
        print(f"  {sample:10} -> {rep}  kept={kept:>6}  dropped_scaffold={dropped}")
        out.append(dst)
    return out


def main():
    args = parse_args()
    outdir = args.outdir
    workdir = outdir / "work"
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] normalising chromosome names ({args.raw_dir})")
    reps = normalise(args.raw_dir, workdir, args.keep_scaffolds)

    print("[2/4] merging strand-aware across replicates")
    allbed = workdir / "all.bed"
    with open(allbed, "w") as fout:
        subprocess.run(["sort", "-k1,1", "-k2,2n"] + [str(p) for p in reps], stdout=fout, check=True)
    merged = workdir / "merged_raw.bed"
    with open(merged, "w") as fout:
        subprocess.run(["bedtools", "merge", "-i", str(allbed), "-s",
                        "-c", "4,5,6", "-o", "distinct,sum,distinct"], stdout=fout, check=True)

    print("[3/4] applying reproducibility filter and collapsing to 1 nt anchors")
    support_hist = Counter()
    widths = []
    rows = []
    with open(merged) as fin:
        for line in fin:
            c = line.rstrip("\n").split("\t")
            chrom, start, end = c[0], int(c[1]), int(c[2])
            reps_here = sorted(set(c[3].split(",")))
            strand = c[5].split(",")[0]
            n = len(reps_here)
            support_hist[n] += 1
            if n < args.min_reps:
                continue
            widths.append(end - start)
            anchor = start if args.anchor == "start" else (start + end) // 2
            rows.append((chrom, anchor, strand, n))

    rows.sort(key=lambda r: (r[0], r[1]))
    # Retain one record per (chromosome, anchor, strand).
    seen = set()
    final = []
    for chrom, anchor, strand, n in rows:
        key = (chrom, anchor, strand)
        if key in seen:
            continue
        seen.add(key)
        final.append((chrom, anchor, strand, n))

    tag = f"min{args.min_reps}rep"
    out_bed = outdir / f"THRAP3_merged_{tag}_anchors.bed"
    with open(out_bed, "w") as fout:
        for i, (chrom, anchor, strand, n) in enumerate(final, 1):
            fout.write("\t".join([chrom, str(anchor), str(anchor + 1),
                                  f"THRAP3_{n}reps_{i}", str(n), strand]) + "\n")

    print("[4/4] summary")
    total = sum(support_hist.values())
    print(f"  merged regions (>=1 rep): {total}")
    for k in sorted(support_hist):
        print(f"    present in {k}/4 reps: {support_hist[k]:>7}")
    cum = sum(v for k, v in support_hist.items() if k >= args.min_reps)
    print(f"  passing >={args.min_reps} reps: {cum} ({100*cum/total:.1f}%)")
    if widths:
        print(f"  merged width: mean={sum(widths)/len(widths):.1f} max={max(widths)}")
    print(f"  anchors written: {len(final)} (deduped from {len(rows)})")
    print(f"  -> {out_bed}")


if __name__ == "__main__":
    main()
