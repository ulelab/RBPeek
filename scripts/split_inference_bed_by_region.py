#!/usr/bin/env python3
"""
Partition an inference BED into exonic and intronic loci, and write the corresponding
normalisation regions for intersect_inference_bed.py (--norm-bed).

Classification is strand-aware and gives exons priority:

  exonic      the anchor lies within an annotated exon on the same strand
  intronic    the anchor lies within a gene on the same strand but in no exon
  intergenic  neither; counted in the summary but not written

The exons of all transcripts are merged before classification, so an anchor that is exonic in
any transcript is classified as exonic, and "intronic" means within a gene and outside every
merged exon. The two sets are therefore disjoint. An anchor within a gene on the opposite
strand is intergenic with respect to that gene.

Outputs
  <prefix>_exonic.bed, <prefix>_intronic.bed   the two locus sets
  regions_exonic.bed                           merged exons, per strand
  regions_intronic.bed                         gene bodies minus merged exons, per strand

--drop-chrM removes mitochondrial anchors before classification. Mitochondrial rRNA is a major
source of background in eCLIP libraries, and intersect_inference_bed.py excludes mitochondrial
peaks from its normalisation.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

CHRM = {"chrM", "chrMT", "MT", "M"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--bed", required=True, type=Path, help="Inference BED (BED6, 1 nt anchors)")
    p.add_argument("-g", "--gtf", required=True, type=Path,
                   help="Gene annotation GTF of the same assembly as the BED (the THRAP3 analysis used "
                        "GENCODE v39, GRCh38)")
    p.add_argument("-o", "--outdir", type=Path, default=Path("THRAP3"), help="Output directory")
    p.add_argument("--prefix", default=None, help="Output basename (default: the input BED's stem)")
    p.add_argument("--drop-chrM", dest="drop_chrm", action="store_true",
                   help="Remove mitochondrial anchors before classifying")
    p.add_argument(
        "--gene-feature",
        default="gene",
        help="GTF feature naming a whole gene body (default 'gene'; use 'transcript' if the GTF has no gene rows)",
    )
    return p.parse_args()


def run(cmd, stdout=None):
    # universal_newlines is used in place of text= for compatibility with Python 3.6.
    r = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, universal_newlines=True)
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
    for tmp in (raw, srt):
        if tmp.exists():
            tmp.unlink()
    return n


def subset(anchors: Path, regions: Path, out: Path, invert: bool = False) -> int:
    flag = "-v" if invert else "-u"
    with open(out, "w") as fh:
        run(["bedtools", "intersect", "-a", str(anchors), "-b", str(regions), "-s", flag], stdout=fh)
    return sum(1 for _ in open(out))


def require(tool):
    if shutil.which(tool) is None:
        sys.exit(f"{tool} not found on PATH")


def main():
    args = parse_args()
    require("bedtools")
    require("sort")
    outdir = args.outdir
    work = outdir / "region_work"
    work.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.bed.stem

    # Match the chromosome naming of the GTF to that of the inference BED.
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

    # Normalisation regions for intersect_inference_bed.py --norm-bed.
    reg_exonic = outdir / "regions_exonic.bed"
    reg_intronic = outdir / "regions_intronic.bed"
    shutil.copyfile(str(exons), str(reg_exonic))
    with open(reg_intronic, "w") as fh:
        run(["bedtools", "subtract", "-s", "-a", str(genes), "-b", str(exons)], stdout=fh)

    anchors = args.bed
    n_input = sum(1 for _ in open(args.bed))
    n_chrm = 0
    if args.drop_chrm:
        anchors = work / "anchors_nochrM.bed"
        with open(args.bed) as fin, open(anchors, "w") as fout:
            for line in fin:
                if line.split("\t", 1)[0] in CHRM:
                    n_chrm += 1
                    continue
                fout.write(line)

    print("[2/3] classifying anchors (strand-aware, exon-priority)")
    total = n_input - n_chrm
    exonic = outdir / f"{prefix}_exonic.bed"
    n_exonic = subset(anchors, exons, exonic)

    genic = work / "genic.bed"
    subset(anchors, genes, genic)
    intronic = outdir / f"{prefix}_intronic.bed"
    n_intronic = subset(genic, exons, intronic, invert=True)

    n_intergenic = total - n_exonic - n_intronic

    print("[3/3] summary")
    print(f"      input anchors      {n_input:>8,}")
    if args.drop_chrm:
        print(f"      chrM dropped       {n_chrm:>8,}")
    print(f"      classified         {total:>8,}")
    print(f"      exonic             {n_exonic:>8,}  ({100*n_exonic/total:5.1f}%)  -> {exonic}")
    print(f"      intronic           {n_intronic:>8,}  ({100*n_intronic/total:5.1f}%)  -> {intronic}")
    print(f"      intergenic         {n_intergenic:>8,}  ({100*n_intergenic/total:5.1f}%)")
    print(f"      normalisation regions -> {reg_exonic}, {reg_intronic}")
    if n_exonic == 0 or n_intronic == 0:
        sys.exit("one subset is empty - check that the GTF assembly matches the inference BED")


if __name__ == "__main__":
    main()
