#!/usr/bin/env python3
"""
Split an inference BED into exonic and intronic subsets so the same loci can be run
through intersect_inference_bed.py twice and the two enrichment profiles compared.

Classification is strand-aware and exon-priority:

  exonic     the anchor falls inside an annotated exon on its own strand
  intronic   the anchor falls inside a gene on its own strand but in no exon
  intergenic neither (written out only with --keep-intergenic, never analysed by default)

Exon-priority matters because an anchor can be exonic in one transcript and intronic in
another. Merging every transcript's exons before classifying resolves that consistently: an
anchor exonic in ANY transcript is exonic here, and intronic means 'inside a gene and in no
merged exon'. The two sets are therefore disjoint by construction and can be compared
without double-counting.

Strandedness is not optional here: an anchor sitting inside a gene on the opposite strand is
intergenic with respect to that gene, and treating it otherwise would put antisense loci in
the intronic set.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--bed", required=True, type=Path, help="Inference BED (BED6, 1 nt anchors)")
    p.add_argument("-g", "--gtf", required=True, type=Path, help="Gene annotation GTF (must match the BED's assembly)")
    p.add_argument("-o", "--outdir", type=Path, default=Path("THRAP3"), help="Output directory")
    p.add_argument("--prefix", default=None, help="Output basename (default: the input BED's stem)")
    p.add_argument("--keep-intergenic", action="store_true", help="Also write the intergenic subset")
    p.add_argument(
        "--gene-feature",
        default="gene",
        help="GTF feature naming a whole gene body (default 'gene'; use 'transcript' if the GTF has no gene rows)",
    )
    return p.parse_args()


def run(cmd, stdout=None):
    r = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(str(c) for c in cmd)}\n{r.stderr[:600]}")
    return r


def gtf_to_bed(gtf: Path, feature: str, out: Path, chr_prefix: bool) -> int:
    """Extract one GTF feature type to sorted, strand-aware-merged BED6."""
    raw = out.with_suffix(".raw")
    add = '($1 ~ /^chr/) ? $1 : "chr"$1' if chr_prefix else "$1"
    awk = (
        f'BEGIN{{FS=OFS="\\t"}} !/^#/ && $3=="{feature}" '
        f'{{c = {add}; if (c=="chrMT") c="chrM"; print c, $4-1, $5, ".", ".", $7}}'
    )
    with open(raw, "w") as fh:
        run(["awk", awk, str(gtf)], stdout=fh)
    n = sum(1 for _ in open(raw))
    if n == 0:
        sys.exit(f"no '{feature}' rows found in {gtf} - check --gene-feature and the GTF's feature column")
    srt = out.with_suffix(".sorted")
    with open(srt, "w") as fh:
        run(["sort", "-k1,1", "-k2,2n", str(raw)], stdout=fh)
    with open(out, "w") as fh:
        run(["bedtools", "merge", "-i", str(srt), "-s", "-c", "6", "-o", "distinct"], stdout=fh)
    # bedtools merge -s emits chrom,start,end,strand; pad back to BED6 with strand in col6
    padded = out.with_suffix(".bed6")
    with open(out) as fin, open(padded, "w") as fout:
        for line in fin:
            c = line.rstrip("\n").split("\t")
            fout.write("\t".join([c[0], c[1], c[2], ".", ".", c[3]]) + "\n")
    padded.replace(out)
    raw.unlink(missing_ok=True)
    srt.unlink(missing_ok=True)
    return n


def subset(anchors: Path, regions: Path, out: Path, invert: bool = False) -> int:
    flag = "-v" if invert else "-u"
    with open(out, "w") as fh:
        run(["bedtools", "intersect", "-a", str(anchors), "-b", str(regions), "-s", flag], stdout=fh)
    return sum(1 for _ in open(out))


def main():
    args = parse_args()
    outdir = args.outdir
    work = outdir / "region_work"
    work.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.bed.stem

    # The inference BED is UCSC-style ("chr1"); Ensembl GTFs are not. Detect and match.
    bed_chr = open(args.bed).readline().split("\t")[0].startswith("chr")
    gtf_chr = False
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            gtf_chr = line.split("\t")[0].startswith("chr")
            break
    need_prefix = bed_chr and not gtf_chr
    print(f"[1/3] reading annotation ({args.gtf})")
    print(f"      inference BED is {'chr-prefixed' if bed_chr else 'Ensembl-style'}; "
          f"GTF is {'chr-prefixed' if gtf_chr else 'Ensembl-style'}"
          f"{' -> adding chr prefix to GTF' if need_prefix else ''}")
    if gtf_chr and not bed_chr:
        sys.exit("GTF is chr-prefixed but the inference BED is not; normalise the BED first.")

    exons = work / "exons.bed"
    genes = work / "genes.bed"
    n_ex = gtf_to_bed(args.gtf, "exon", exons, need_prefix)
    n_gn = gtf_to_bed(args.gtf, args.gene_feature, genes, need_prefix)
    print(f"      {n_ex:,} exon rows -> {sum(1 for _ in open(exons)):,} merged intervals")
    print(f"      {n_gn:,} {args.gene_feature} rows -> {sum(1 for _ in open(genes)):,} merged intervals")

    print("[2/3] classifying anchors (strand-aware, exon-priority)")
    total = sum(1 for _ in open(args.bed))
    exonic = outdir / f"{prefix}_exonic.bed"
    n_exonic = subset(args.bed, exons, exonic)

    genic = work / "genic.bed"
    subset(args.bed, genes, genic)
    intronic = outdir / f"{prefix}_intronic.bed"
    n_intronic = subset(genic, exons, intronic, invert=True)

    n_intergenic = total - n_exonic - n_intronic
    if args.keep_intergenic:
        intergenic = outdir / f"{prefix}_intergenic.bed"
        subset(args.bed, genes, intergenic, invert=True)
        print(f"      intergenic -> {intergenic}")

    print("[3/3] summary")
    print(f"      total anchors      {total:>8,}")
    print(f"      exonic             {n_exonic:>8,}  ({100*n_exonic/total:5.1f}%)  -> {exonic}")
    print(f"      intronic           {n_intronic:>8,}  ({100*n_intronic/total:5.1f}%)  -> {intronic}")
    print(f"      intergenic         {n_intergenic:>8,}  ({100*n_intergenic/total:5.1f}%)")
    print("      exon-priority: all transcripts' exons are merged first, so an anchor that is")
    print("      exonic in ANY transcript is exonic here. The two sets are disjoint by")
    print("      construction - intronic is 'genic AND in no merged exon'.")
    if n_exonic == 0 or n_intronic == 0:
        sys.exit("one subset is empty - check that the GTF assembly matches the inference BED")


if __name__ == "__main__":
    main()
