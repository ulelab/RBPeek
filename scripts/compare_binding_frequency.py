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
    p.add_argument("-o", "--out", default=None, help="Write the full table here as TSV")
    p.add_argument("--top", type=int, default=25, help="How many rows to print (default 25)")
    return p.parse_args()


def load(path):
    """Return (n_loci, {protein: n_loci_with_signal}) from a binf_summary.tsv."""
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
        counts = dict((name, 0) for _, name in idx)
        n = 0
        for row in reader:
            if not row:
                continue
            n += 1
            for i, name in idx:
                # float() then compare, since these are written in %g form
                if row[i] not in ("0", "0.0", "") and float(row[i]) > 0:
                    counts[name] += 1
    return n, counts


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
    nA, cA = load(args.a)
    nB, cB = load(args.b)
    shared = sorted(set(cA) & set(cB))
    only_a = sorted(set(cA) - set(cB))
    only_b = sorted(set(cB) - set(cA))
    print("%s: %d loci, %d panel columns" % (args.label_a, nA, len(cA)))
    print("%s: %d loci, %d panel columns" % (args.label_b, nB, len(cB)))
    if only_a or only_b:
        print("WARNING: panel differs between the two runs; comparing the %d shared columns only"
              % len(shared))
        for name, lst in ((args.label_a, only_a), (args.label_b, only_b)):
            if lst:
                print("  only in %s: %s%s" % (name, ", ".join(lst[:5]), " ..." if len(lst) > 5 else ""))

    rows = []
    for p in shared:
        xa, xb = cA[p], cB[p]
        if max(xa, xb) < args.min_loci:
            continue
        fa, fb = xa / float(nA), xb / float(nB)
        rows.append((p, fa, fb, fa - fb, z_two_proportion(xa, nA, xb, nB), xa, xb))
    if not rows:
        sys.exit("no protein cleared --min-loci %d in either set" % args.min_loci)
    rows.sort(key=lambda r: -r[4])

    hdr = ("protein", "frac_" + args.label_a, "frac_" + args.label_b, "diff", "z",
           "n_" + args.label_a, "n_" + args.label_b)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for r in rows:
                fh.write("%s\t%.4f\t%.4f\t%+.4f\t%+.2f\t%d\t%d\n" % r)
        print("\nwrote %d rows -> %s" % (len(rows), args.out))

    print("\n=== most enriched in %s relative to %s ===" % (args.label_a, args.label_b))
    print("%-30s %9s %9s %9s %8s" % ("protein", "frac_A", "frac_B", "diff", "z"))
    for r in rows[: args.top]:
        print("%-30s %9.3f %9.3f %+9.3f %+8.1f" % (r[0], r[1], r[2], r[3], r[4]))
    print("\n=== most enriched in %s (i.e. background-associated) ===" % args.label_b)
    for r in rows[-5:]:
        print("%-30s %9.3f %9.3f %+9.3f %+8.1f" % (r[0], r[1], r[2], r[3], r[4]))
    print("\nz is a two-proportion test on locus counts. It scales with set size, so treat it")
    print("as a ranking aid, not a p-value: with thousands of loci, small differences reach")
    print("large z. The effect size that matters is `diff`.")


if __name__ == "__main__":
    main()
