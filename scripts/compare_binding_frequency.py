#!/usr/bin/env python3
"""
Compare how often each panel protein binds two different inference-locus sets.

Written for proximity-labelling data, where the positional metrics RBPeek normally ranks on
carry little protein-specific information: biotin labels within a spatial radius, so a peak's
position reflects distance from the bait compartment rather than any RBP's own footprint.
What does survive is WHICH loci a protein binds, so this compares binding frequency instead.

For each protein present in both tables:

    frac_A = fraction of set-A loci where the protein has any signal
    frac_B = same for set B
    diff   = frac_A - frac_B

with a two-proportion z-test. A protein binding both sets equally is reporting background;
a protein binding A much more often is a candidate specific to whatever defines A.

Typical use, after running both halves of run_centro_controlled.sbatch:

    python3 scripts/compare_binding_frequency.py \\
      -a results/centro_specific/binf_summary.tsv   --label-a centro_specific \\
      -b results/centro_ntcontrol/binf_summary.tsv  --label-b NT_control \\
      -o results/centro_specific/binding_frequency_vs_control.tsv

Python 3.6 compatible so it runs with the login node's system python.
"""

import argparse
import csv
import math
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-a", required=True, help="binf_summary.tsv for set A (the set of interest)")
    p.add_argument("-b", required=True, help="binf_summary.tsv for set B (the control set)")
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--min-loci", type=int, default=20,
                   help="Require this many signal-bearing loci in at least one set (default 20)")
    p.add_argument("--match-support", action="store_true",
                   help="Decile-match the two locus sets on total panel support before comparing. "
                        "Strongly recommended whenever one set was DERIVED from the other (e.g. by "
                        "subtracting controls), because that changes the loci's bindability and "
                        "confounds a raw frequency comparison.")
    p.add_argument("--effect", choices=["h", "diff"], default="h",
                   help="Ranking metric. 'h' (default) is Cohen's h, a scale-free effect size for "
                        "proportions; 'diff' is the raw difference, which favours high-coverage "
                        "columns because a protein binding 5%% of loci cannot show a large one.")
    p.add_argument("--seed", type=int, default=0, help="Seed for match sampling (default 0)")
    p.add_argument("-o", "--out", default=None, help="Write the full table here as TSV")
    p.add_argument("--top", type=int, default=25, help="How many rows to print (default 25)")
    return p.parse_args()


def load(path):
    """Return (names, rows, support) where rows[i] is a per-protein 0/1 presence list and
    support[i] is that locus's summed signal across the whole panel."""
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.exit("empty file: %s" % path)
        idx = [(i, c[: -len("_total_overlaps")]) for i, c in enumerate(header)
               if c.endswith("_total_overlaps")]
        if not idx:
            sys.exit("no *_total_overlaps columns in %s" % path)
        names = [name for _, name in idx]
        rows, support = [], []
        for row in reader:
            if not row:
                continue
            pres, tot = [], 0.0
            for i, _ in idx:
                v = 0.0
                if row[i] not in ("0", "0.0", ""):
                    v = float(row[i])
                pres.append(1 if v > 0 else 0)
                tot += v
            rows.append(pres)
            support.append(tot)
    return names, rows, support


def match_on_support(supA, supB, seed, nbins=10):
    """
    Return (idxA, idxB): indices of loci sampled so both sets share a support distribution.

    Necessary when one set is derived from the other. Subtracting control peaks from a bait
    set removes the strongly-bound, accessible loci, so the remainder is a systematically
    lower-support population and EVERY protein then looks control-enriched. On the centrosome
    data the mean binding frequency was 0.373 against 0.561 - a 0.665 global ratio that no
    per-protein statistic can see past, and which pushes sparse columns to the top purely
    because they cannot go far negative.
    """
    import random
    rng = random.Random(seed)
    pooled = sorted(supA + supB)
    edges = [pooled[int(round(q * (len(pooled) - 1) / float(nbins)))] for q in range(1, nbins)]

    def binof(v):
        lo = 0
        for e in edges:
            if v <= e:
                return lo
            lo += 1
        return lo

    binsA, binsB = [binof(v) for v in supA], [binof(v) for v in supB]
    keepA, keepB = [], []
    for b in range(nbins):
        ia = [i for i, x in enumerate(binsA) if x == b]
        ib = [i for i, x in enumerate(binsB) if x == b]
        take = min(len(ia), len(ib))
        if not take:
            continue
        keepA.extend(ia if len(ia) == take else rng.sample(ia, take))
        keepB.extend(ib if len(ib) == take else rng.sample(ib, take))
    return keepA, keepB


def cohens_h(p1, p2):
    """Scale-free effect size for two proportions; behaves sensibly near 0 and 1."""
    p1 = min(max(p1, 0.0), 1.0)
    p2 = min(max(p2, 0.0), 1.0)
    return 2 * math.asin(math.sqrt(p1)) - 2 * math.asin(math.sqrt(p2))


def z_two_proportion(x1, n1, x2, n2):
    """Two-proportion z-test; returns 0.0 when undefined rather than blowing up."""
    if n1 == 0 or n2 == 0:
        return 0.0
    p1, p2 = x1 / float(n1), x2 / float(n2)
    p = (x1 + x2) / float(n1 + n2)
    se = math.sqrt(p * (1 - p) * (1.0 / n1 + 1.0 / n2))
    if se == 0:
        return 0.0
    return (p1 - p2) / se


def main():
    args = parse_args()
    namesA, rowsA, supA = load(args.a)
    namesB, rowsB, supB = load(args.b)
    nA, nB = len(rowsA), len(rowsB)
    print("%s: %d loci, %d panel columns" % (args.label_a, nA, len(namesA)))
    print("%s: %d loci, %d panel columns" % (args.label_b, nB, len(namesB)))

    posA = dict((n, i) for i, n in enumerate(namesA))
    posB = dict((n, i) for i, n in enumerate(namesB))
    shared = sorted(set(namesA) & set(namesB))
    only_a = sorted(set(namesA) - set(namesB))
    only_b = sorted(set(namesB) - set(namesA))
    if only_a or only_b:
        print("WARNING: panel differs between the two runs; comparing the %d shared columns only"
              % len(shared))
        for name, lst in ((args.label_a, only_a), (args.label_b, only_b)):
            if lst:
                print("  only in %s: %s%s" % (name, ", ".join(lst[:5]), " ..." if len(lst) > 5 else ""))

    def freqs(rows, keep, pos):
        m = len(keep) if keep is not None else len(rows)
        out = {}
        for p in shared:
            j = pos[p]
            if keep is None:
                c = sum(r[j] for r in rows)
            else:
                c = sum(rows[i][j] for i in keep)
            out[p] = (c, m)
        return out

    def global_ratio(fa, fb):
        ma = sum(c / float(m) for c, m in fa.values()) / len(fa)
        mb = sum(c / float(m) for c, m in fb.values()) / len(fb)
        return ma, mb

    keepA = keepB = None
    fa, fb = freqs(rowsA, None, posA), freqs(rowsB, None, posB)
    ma, mb = global_ratio(fa, fb)
    print("\nmean binding frequency across the shared panel: %s=%.3f  %s=%.3f  ratio=%.3f"
          % (args.label_a, ma, args.label_b, mb, ma / mb if mb else float("nan")))
    if abs(math.log(ma / mb)) > 0.15 and not args.match_support:
        print("  NOTE this is a large global offset. Every protein will look enriched in the")
        print("  denser set, and sparse columns will float to the top simply because they")
        print("  cannot go far in the other direction. Consider --match-support.")

    if args.match_support:
        keepA, keepB = match_on_support(supA, supB, args.seed)
        fa, fb = freqs(rowsA, keepA, posA), freqs(rowsB, keepB, posB)
        ma2, mb2 = global_ratio(fa, fb)
        print("support-matched: %d vs %d loci; mean frequency now %.3f vs %.3f (ratio %.3f)"
              % (len(keepA), len(keepB), ma2, mb2, ma2 / mb2 if mb2 else float("nan")))
        nA, nB = len(keepA), len(keepB)

    rows = []
    for p in shared:
        xa, ta = fa[p]
        xb, tb = fb[p]
        if max(xa, xb) < args.min_loci:
            continue
        pa, pb = xa / float(ta), xb / float(tb)
        rows.append((p, pa, pb, pa - pb, cohens_h(pa, pb), z_two_proportion(xa, ta, xb, tb), xa, xb))
    if not rows:
        sys.exit("no protein cleared --min-loci %d in either set" % args.min_loci)
    key = 4 if args.effect == "h" else 3
    rows.sort(key=lambda r: -r[key])
    metric = "cohens_h" if args.effect == "h" else "diff"

    hdr = ("protein", "frac_" + args.label_a, "frac_" + args.label_b, "diff", "cohens_h", "z",
           "n_" + args.label_a, "n_" + args.label_b)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for r in rows:
                fh.write("%s\t%.4f\t%.4f\t%+.4f\t%+.4f\t%+.2f\t%d\t%d\n" % r)
        print("\nwrote %d rows -> %s (sorted by %s)" % (len(rows), args.out, metric))

    print("\n=== most enriched in %s relative to %s (by %s) ==="
          % (args.label_a, args.label_b, metric))
    print("%-30s %9s %9s %9s %9s" % ("protein", "frac_A", "frac_B", "diff", "cohen_h"))
    for r in rows[: args.top]:
        print("%-30s %9.3f %9.3f %+9.3f %+9.3f" % (r[0], r[1], r[2], r[3], r[4]))
    print("\n=== most enriched in %s ===" % args.label_b)
    for r in rows[-5:]:
        print("%-30s %9.3f %9.3f %+9.3f %+9.3f" % (r[0], r[1], r[2], r[3], r[4]))
    print("\nz is a two-proportion test and scales with set size, so treat it as a ranking aid")
    print("rather than a p-value. cohens_h is the scale-free effect size; `diff` favours")
    print("high-coverage columns, since a protein binding 5% of loci cannot show a large one.")


if __name__ == "__main__":
    main()
