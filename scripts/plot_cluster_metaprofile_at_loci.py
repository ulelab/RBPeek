#!/usr/bin/env python3
"""
Metaprofile of one heatmap cluster's proteins, computed around a DIFFERENT set of loci.

The motivating question: cluster 1 of the centro-specific run is defined by a particular set
of RBPs - do those same RBPs show the same binding profile around the CONTROL peaks?

Why this needs its own run rather than a read of existing outputs
----------------------------------------------------------------
A metaprofile is mean support at every offset in -window..+window. That vector lives in the
(n_loci x 2*window+1) counts matrix that intersect_inference_bed.py builds in memory and
never writes out. binf_summary.tsv keeps only five SUMMARY statistics per protein per locus
(total_overlaps, variance, skew, kurtosis, max_binding_offset) - how MUCH a protein binds
each locus, not WHERE inside the window. So the profile cannot be recovered from it; the
counts have to be recomputed against the target loci.

This script recomputes them for a handful of proteins only, so it takes seconds rather than
re-running the whole panel. It imports the arithmetic from intersect_inference_bed.py rather
than reimplementing it, so slop, strand-aware intersection, the strand flip and --panel-anchor
behave identically to a full run.

Choosing "cluster 1 proteins"
-----------------------------
binf_heatmap_clusters.tsv assigns LOCI to clusters; it has no protein column. A cluster's
proteins are therefore derived: those most enriched in presence at that cluster's loci
relative to the other clusters. Pass --proteins to override with an explicit list.

Needs the rbpeek conda env (bedtools, and the panel BEDs under --xldir).
"""

import argparse
import csv
import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Categorical slots 1-8 of the reference theme, in their validated order.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a63d2", "#00a3b5", "#c2477f", "#6b7280"]
INK, MUTED, GRID = "#1a1a19", "#52514e", "#e6e5e0"


def load_engine():
    here = Path(__file__).resolve().parent / "intersect_inference_bed.py"
    spec = importlib.util.spec_from_file_location("iib", str(here))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True, help="binf_heatmap_clusters.tsv defining the cluster")
    p.add_argument("--cluster-summary", required=True,
                   help="binf_summary.tsv from the SAME run, used to work out which proteins mark the cluster")
    p.add_argument("--cluster", default="1")
    p.add_argument("--n-proteins", type=int, default=6)
    p.add_argument("--proteins", default=None, help="Comma-separated names, overriding the automatic pick")
    p.add_argument("-b", "--bed", required=True, help="Loci to profile AROUND (e.g. the control anchors BED)")
    p.add_argument("-x", "--xldir", required=True, help="Panel root, as passed to intersect_inference_bed.py")
    p.add_argument("-s", "--samplesheet", required=True, help="Panel samplesheet (file, group)")
    p.add_argument("--window", type=int, default=100)
    p.add_argument("--panel-anchor", choices=["start", "midpoint"], default="midpoint")
    p.add_argument("--gaussian-sigma", type=float, default=2.0)
    p.add_argument("--genome", default=None, help="Genome sizes file (defaults to the engine's DEFAULT_GENOME)")
    p.add_argument("--label", default="control loci", help="What the target loci are, for the title")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--out-table", default=None)
    return p.parse_args()


def cluster_proteins(clusters_tsv, summary_tsv, cluster, n, min_loci=30):
    """Proteins most enriched in presence at the cluster's loci versus the other clusters."""
    lab = {}
    for r in csv.DictReader(open(clusters_tsv), delimiter="\t"):
        lab[r["binf_chr_start_end"]] = r["heatmap_cluster"]

    hdr = open(summary_tsv).readline().rstrip("\n").split("\t")
    keyi = hdr.index("binf_chr_start_end")
    cols = [(i, c[: -len("_total_overlaps")]) for i, c in enumerate(hdr) if c.endswith("_total_overlaps")]
    names = [n_ for _, n_ in cols]
    inC, other = [], []
    with open(summary_tsv) as fh:
        next(fh)
        for line in fh:
            r = line.rstrip("\n").split("\t")
            L = lab.get(r[keyi], "NA")
            pres = [0 if r[i] in ("0", "0.0", "") or float(r[i]) <= 0 else 1 for i, _ in cols]
            if L == cluster:
                inC.append(pres)
            elif L != "NA":
                other.append(pres)
    if not inC:
        sys.exit("cluster %r not present in %s" % (cluster, clusters_tsv))
    A, B = np.array(inC), np.array(other)
    print("cluster %s: %d loci; other clusters: %d loci" % (cluster, len(A), len(B)))
    scored = []
    for j, nm in enumerate(names):
        a = A[:, j].mean()
        b = B[:, j].mean() if len(B) else 0.0
        if b > 0 and A[:, j].sum() >= min_loci:
            scored.append((a / b, nm, a, b))
    scored.sort(reverse=True)
    print("\nproteins marking cluster %s (presence at cluster loci vs other clusters):" % cluster)
    for r, nm, a, b in scored[:n]:
        print("  %-28s %.3f vs %.3f   ratio %.2f" % (nm, a, b, r))
    return [s[1] for s in scored[:n]]


def main():
    args = parse_args()
    eng = load_engine()
    genome = args.genome or eng.DEFAULT_GENOME
    if not Path(genome).exists():
        sys.exit("genome sizes file not found: %s (pass --genome)" % genome)

    picked = ([p.strip() for p in args.proteins.split(",")] if args.proteins
              else cluster_proteins(args.clusters, args.cluster_summary, args.cluster, args.n_proteins))

    sources = eng.uniquify_names(eng.load_samplesheet_inputs(Path(args.samplesheet), Path(args.xldir)))
    by_name = dict(sources)
    missing = [p for p in picked if p not in by_name]
    if missing:
        sys.exit("not in the samplesheet: %s" % ", ".join(missing))

    binf = Path(args.bed)
    print("\nprofiling %d protein(s) around %s (%d loci)"
          % (len(picked), binf.name, sum(1 for _ in open(str(binf)))))

    with tempfile.TemporaryDirectory(prefix="submeta_") as tmp:
        tmpdir = Path(tmp)
        keys, index, windows = eng.load_binf_and_prepare_windows(
            binf_path=binf, window=args.window, genome=genome, tmpdir=tmpdir)
        n_binf = len(keys)
        # same chromosome-convention guard the main tool applies
        fixed = eng.harmonise_panel_chroms([(p, by_name[p]) for p in picked], binf, tmpdir)
        offsets = np.arange(-args.window, args.window + 1)
        profiles = {}
        for name, path in fixed:
            counts = eng.compute_counts_for_protein(
                merged_xl_bed=path, binf_windows_bed=windows, binf_index=index,
                window=args.window, n_binf=n_binf, panel_anchor=args.panel_anchor)
            profiles[name] = eng.smooth_metaprofile_gaussian(
                counts.mean(axis=0), sigma=args.gaussian_sigma)
            print("  %-28s mean support at offset 0 = %.3f" % (name, profiles[name][args.window]))

    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    for i, name in enumerate(picked):
        ax.plot(offsets, profiles[name], color=PALETTE[i % len(PALETTE)], linewidth=2,
                label=name, solid_capstyle="round")
    ax.axvline(0, color=MUTED, linewidth=1, alpha=0.35, zorder=0)
    ax.set_xlim(-args.window, args.window)
    ax.set_xlabel("Relative nucleotide position (nt)", color=MUTED)
    ax.set_ylabel("Mean crosslink support per locus", color=MUTED)
    ax.set_title("Cluster-%s proteins profiled around %s" % (args.cluster, args.label), color=INK)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=MUTED)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5),
              labelcolor=MUTED, fontsize=8)
    fig.tight_layout(rect=[0, 0, 0.99, 1])
    d = os.path.dirname(args.out)
    if d:
        os.makedirs(d, exist_ok=True)
    fig.savefig(args.out, dpi=200, facecolor="#fcfcfb", bbox_inches="tight")
    print("\nwrote %s" % args.out)

    if args.out_table:
        with open(args.out_table, "w") as fh:
            fh.write("offset\t" + "\t".join(picked) + "\n")
            for k, off in enumerate(offsets):
                fh.write(str(off) + "\t" + "\t".join("%.6g" % profiles[p][k] for p in picked) + "\n")
        print("wrote %s" % args.out_table)


if __name__ == "__main__":
    main()
