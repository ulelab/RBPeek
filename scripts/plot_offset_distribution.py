#!/usr/bin/env python3
"""
Where do a set of proteins bind, relative to two different sets of inference loci?

Reads two binf_summary.tsv files and plots, per protein, the distribution of
`max_binding_offset` across loci in each set. Proteins are chosen by how enriched they are
in a given heatmap row cluster of the first set.

IMPORTANT - this is NOT a metaprofile. binf_summary.tsv stores per-locus SUMMARY statistics,
not the per-nucleotide count vectors, so the true +/-window profile cannot be reconstructed
from it. What is plotted is the distribution of each locus's single strongest offset. A
metaprofile would need intersect_inference_bed.py re-run with these loci as the inference BED.

Loci where a protein has no signal are excluded. compute_summary_stats takes an argmax over
an all-zero vector for those, which returns offset 0 - so including them piles a large fake
spike onto the centre (verified: 100% of zero-signal loci report offset 0).
"""

import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Slots 1 and 2 of the reference categorical theme. Validated rather than eyeballed:
# worst-of-protan/deutan OKLab dE = 24.7 (target >=8), normal-vision dE = 33.6 (floor 15),
# contrast vs the #fcfcfb surface 4.30:1 and 3.12:1.
COLOR_A, COLOR_B = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#1a1a19", "#52514e", "#e6e5e0"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-a", required=True, help="binf_summary.tsv for set A (the cluster source)")
    p.add_argument("-b", required=True, help="binf_summary.tsv for set B (e.g. the control loci)")
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--clusters", required=True, help="binf_heatmap_clusters.tsv for set A")
    p.add_argument("--cluster", default="1", help="Which row cluster to take proteins from (default 1)")
    p.add_argument("--n-proteins", type=int, default=9, help="How many cluster-enriched proteins to plot")
    p.add_argument("--proteins", default=None, help="Comma-separated protein names, overriding the automatic pick")
    p.add_argument("--window", type=int, default=100, help="Plot +/- this many nt (default 100, the full --window of the run)")
    p.add_argument("--bin", type=int, default=10, help="Offset bin width in nt (default 10)")
    p.add_argument("--min-loci", type=int, default=30, help="Skip a protein with fewer signal-bearing loci")
    p.add_argument("-o", "--out", required=True)
    return p.parse_args()


def load(path):
    hdr = open(path).readline().rstrip("\n").split("\t")
    tot = dict((c[: -len("_total_overlaps")], i) for i, c in enumerate(hdr) if c.endswith("_total_overlaps"))
    off = dict((c[: -len("_max_binding_offset")], i) for i, c in enumerate(hdr) if c.endswith("_max_binding_offset"))
    keyi = hdr.index("binf_chr_start_end")
    keys, rows = [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            r = line.rstrip("\n").split("\t")
            keys.append(r[keyi])
            rows.append(r)
    return tot, off, keys, rows


def series(tot_i, off_i, rows, mask=None):
    """Offsets at loci where this protein has signal (zero-signal loci excluded)."""
    out = []
    for n, r in enumerate(rows):
        if mask is not None and not mask[n]:
            continue
        if r[tot_i] in ("0", "0.0", ""):
            continue
        if float(r[tot_i]) <= 0:
            continue
        out.append(float(r[off_i]))
    return np.array(out)


def main():
    args = parse_args()
    totA, offA, keysA, rowsA = load(args.a)
    totB, offB, keysB, rowsB = load(args.b)

    cl = {}
    for r in csv.DictReader(open(args.clusters), delimiter="\t"):
        cl[r["binf_chr_start_end"]] = r["heatmap_cluster"]
    labs = np.array([cl.get(k, "NA") for k in keysA])
    inC = labs == args.cluster
    other = np.isin(labs, [c for c in set(labs) if c not in ("NA", args.cluster)])
    if inC.sum() == 0:
        sys.exit("cluster %r not found in %s" % (args.cluster, args.clusters))
    print("set A: %d loci, cluster %s = %d, other clusters = %d"
          % (len(keysA), args.cluster, inC.sum(), other.sum()))
    print("set B: %d loci" % len(keysB))

    shared = [p for p in totA if p in totB]
    if args.proteins:
        picked = [p.strip() for p in args.proteins.split(",")]
        missing = [p for p in picked if p not in shared]
        if missing:
            sys.exit("not in both tables: %s" % ", ".join(missing))
    else:
        # enrichment of presence inside the cluster versus the other clusters
        scored = []
        for p in shared:
            ti = totA[p]
            pres = np.array([0.0 if r[ti] in ("0", "0.0", "") else (1.0 if float(r[ti]) > 0 else 0.0)
                             for r in rowsA])
            a, b = pres[inC].mean(), pres[other].mean() if other.sum() else 0.0
            if b > 0 and pres[inC].sum() >= args.min_loci:
                scored.append((a / b, p, a, b))
        scored.sort(reverse=True)
        picked = [s[1] for s in scored[: args.n_proteins]]
        print("\nproteins most enriched in cluster %s (presence ratio vs other clusters):" % args.cluster)
        for r, p, a, b in scored[: args.n_proteins]:
            print("  %-28s %.3f vs %.3f   ratio %.2f" % (p, a, b, r))

    W, BW = args.window, args.bin
    edges = np.arange(-W - BW / 2.0, W + BW, BW)
    centres = (edges[:-1] + edges[1:]) / 2.0

    ncol = 3
    nrow = int(np.ceil(len(picked) / float(ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.5 * nrow),
                             sharex=True, squeeze=False)
    for k, p in enumerate(picked):
        ax = axes[k // ncol][k % ncol]
        sa = series(totA[p], offA[p], rowsA, inC)
        sb = series(totB[p], offB[p], rowsB)
        for s, col, lab in ((sa, COLOR_A, "%s cluster %s" % (args.label_a, args.cluster)),
                            (sb, COLOR_B, args.label_b)):
            if len(s) < args.min_loci:
                continue
            # Restrict rather than clip. np.clip would pile every offset beyond the
            # window into the end bins and manufacture edge spikes that look like signal.
            inw = s[(s >= -W) & (s <= W)]
            if len(inw) < args.min_loci:
                continue
            h, _ = np.histogram(inw, bins=edges)
            ax.plot(centres, 100.0 * h / max(h.sum(), 1), color=col, linewidth=2,
                    label="%s (n=%d, %.0f%% in view)" % (lab, len(s), 100.0 * len(inw) / len(s)),
                    solid_capstyle="round")
        ax.axvline(0, color=MUTED, linewidth=1, alpha=0.35, zorder=0)
        ax.set_title(p, fontsize=8, color=INK)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=7)
        ax.legend(frameon=False, fontsize=6, loc="upper right", labelcolor=MUTED)
    for k in range(len(picked), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    for c in range(ncol):
        axes[nrow - 1][c].set_xlabel("offset of strongest binding (nt)", fontsize=8, color=MUTED)
    for r in range(nrow):
        axes[r][0].set_ylabel("% of bound loci", fontsize=8, color=MUTED)
    fig.suptitle("Where cluster-%s proteins bind, relative to each locus set" % args.cluster,
                 fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=200, facecolor="#fcfcfb")
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
