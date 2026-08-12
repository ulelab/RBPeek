#!/usr/bin/env python3
"""
Build a THRAP3 inference BED (`binf`) for intersect_inference_bed.py from the four
Clippy peak calls produced by the Flow CLIP-Seq run (project HA_THRAP3_CLIP,
788995297969977723, execution 188723932105002059, CLIP-Seq v1.7, GRCh38).

Three things have to happen before these peaks can be used as an inference BED:

1. Chromosome naming. Flow ran Clippy against Homo_sapiens.GRCh38.fasta.fai, so the
   peaks are Ensembl-style ("1"). The RBPeek panel under --xldir and the existing
   decoys BED are UCSC-style ("chr1"). intersect_inference_bed.py only normalises
   names on its internal merge path, never for the inference BED, so an unnormalised
   THRAP3 BED would silently intersect nothing.

2. Replicate reproducibility. The four libraries differ ~5x in depth (9,455 to 49,477
   peaks), so a plain union is dominated by whichever library was sequenced deepest.
   Overlapping peaks are merged strand-aware and kept only when at least
   --min-reps distinct replicates contribute.

3. Single-nucleotide anchors. intersect_inference_bed.py computes
   offset = xl_start - binf_site_start (see intersect_inference_bed.py:356), so every
   inference locus must be 1 bp wide or the metaprofile smears by the locus width.
   Each merged region is collapsed to its midpoint.

Output BED6:
    chrom  start  end(=start+1)  name  score  strand
where name is THRAP3_<n>reps_<i> and score is the number of contributing replicates,
so downstream tables can be split by reproducibility tier.
"""

import argparse
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_RAW = REPO / "THRAP3" / "raw"
# Flow sample name -> replicate label used in the peak filenames.
SAMPLE_TO_REP = {
    "THRAP3_1": "R1",
    "THRAP3_2": "R2",
    "THRAP3_3": "R3",
    "THRAP3_L": "R4",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW,
                   help="Directory holding the four *_genome.*_Peaks.bed files from Flow")
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


def normalise(raw_dir: Path, workdir: Path, keep_scaffolds: bool) -> list[Path]:
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
                # column 4 carries the replicate tag; bedtools merge -o distinct on it
                # is what yields the per-region replicate support count.
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
    # Collapsing to midpoints can make two nearby merged regions land on the same
    # anchor; keep one row per (chrom, anchor, strand) so loci stay unique.
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
