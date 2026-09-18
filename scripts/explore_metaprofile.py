#!/usr/bin/env python3
"""
Try metaprofile and ranking options straight from metaprofile_matrix.npz, without rerunning
intersect_inference_bed.py (write the file with its --save-matrix flag).

  python3 scripts/explore_metaprofile.py results/thrap3_exonic/metaprofile_matrix.npz \\
      --curve log --rank-by area -o results/thrap3_exonic/explore_log_area.pdf

Every option keeps the region normalisation (support x 1e6 / region_cdna). For each sample a
curve is built over the offsets, the samples are ranked, and the top N are drawn with the
engine's own plotting code.

--curve
  log           mean over all loci of log1p(support) at each offset  (the engine's curve)
  linear        mean over all loci of support at each offset
  windowlog     mean over all loci of log1p(sum of support within +/-central-window of each offset)
  windowlinear  mean over all loci of that window sum, no log

--rank-by
  area     area under the unsmoothed curve within +/-central-window  (the engine's score for --curve log)
  height0  the unsmoothed curve's value at offset 0
  file     keep the ranking stored in the matrix file

The printout reports how well peak height follows the chosen ranking (pairs of plotted samples
whose smoothed peak heights are out of rank order), and how concentrated each sample's central
signal is (share coming from its top 1% of loci).
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

ENGINE = Path(__file__).resolve().parent / "intersect_inference_bed.py"
YLABELS = {
    "log": "Mean log1p support per locus\n(per M region cDNA{sm})",
    "linear": "Mean support per locus\n(per M region cDNA{sm})",
    "windowlog": "Mean log1p support in the ±{cw} nt window at each offset\n(per M region cDNA{sm})",
    "windowlinear": "Mean support in the ±{cw} nt window at each offset\n(per M region cDNA{sm})",
}


def load_engine():
    spec = importlib.util.spec_from_file_location("rbpeek_engine", ENGINE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("matrix", type=Path, help="metaprofile_matrix.npz written by --save-matrix")
    p.add_argument("--curve", choices=sorted(YLABELS), default="log")
    p.add_argument("--rank-by", choices=["area", "height0", "file"], default="area")
    p.add_argument("--subtract-flank", action="store_true",
                   help="Subtract each curve's mean beyond +/-central-window (after ranking)")
    p.add_argument("--top", type=int, default=10, help="Number of samples to draw (default 10)")
    p.add_argument("--sigma", type=float, default=2.0, help="Gaussian smoothing for display; 0 = none (default 2)")
    p.add_argument("--central-window", type=int, default=None, help="Override the file's central window (nt)")
    p.add_argument("-o", "--out", type=Path, default=None, help="Output PDF or PNG (default: next to the matrix)")
    return p.parse_args()


def build_curve(counts, scale, kind, half):
    x = counts.astype(np.float64) * scale
    if kind in ("windowlog", "windowlinear"):
        # Sum over +/-half around every offset; windows running past the edge are truncated.
        cs = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(x, axis=1)], axis=1)
        n = x.shape[1]
        lo = np.clip(np.arange(n) - half, 0, n)
        hi = np.clip(np.arange(n) + half + 1, 0, n)
        x = np.clip(cs[:, hi] - cs[:, lo], 0.0, None)
    return (np.log1p(x) if kind in ("log", "windowlog") else x).mean(axis=0)


def main():
    args = parse_args()
    eng = load_engine()
    meta, _ = eng.load_counts_matrix(args.matrix, samples=[])
    names = [str(s) for s in meta["samples"]]
    offsets = meta["offsets"].astype(np.int64)
    window = int(offsets[-1])
    cw = int(args.central_window if args.central_window is not None else meta["central_window"])
    central = np.abs(offsets) <= cw
    zero = int(np.where(offsets == 0)[0][0])

    curves, conc = {}, {}
    for i, pn in enumerate(names):
        reg = float(meta["region_cdna"][i])
        if not reg > 0:
            continue
        _, dense = eng.load_counts_matrix(args.matrix, samples=[pn])
        counts = dense[pn]
        curves[pn] = build_curve(counts, 1e6 / reg, args.curve, cw)
        c = np.sort(counts[:, central].sum(axis=1))[::-1]
        k = max(1, int((c > 0).sum()) // 100)
        conc[pn] = float(c[:k].sum() / c.sum()) if c.sum() > 0 else float("nan")

    old_rank = {pn: int(r) for pn, r in zip(names, meta["rank"])}
    if args.rank_by == "file":
        score = {pn: -old_rank[pn] for pn in curves}
    elif args.rank_by == "height0":
        score = {pn: float(curves[pn][zero]) for pn in curves}
    else:
        score = {pn: float(curves[pn][central].sum()) for pn in curves}
    ranked = sorted(curves, key=lambda pn: score[pn], reverse=True)
    top = ranked[:args.top]

    smooth = (lambda v: eng.smooth_metaprofile_gaussian(v, args.sigma)) if args.sigma > 0 else (lambda v: v)
    shown = {}
    for pn in top:
        v = smooth(curves[pn])
        if args.subtract_flank:
            v = v - v[~central].mean()
        shown[pn] = v
    peaks = [float(shown[pn][central].max()) for pn in top]
    inversions = sum(1 for i in range(len(peaks)) for j in range(i + 1, len(peaks)) if peaks[j] > peaks[i])
    pairs = len(peaks) * (len(peaks) - 1) // 2

    print(f"{args.matrix}: {len(meta['loci']):,} loci, {len(curves)} samples with region cDNA")
    print(f"curve = {args.curve}, rank by = {args.rank_by}, central window = ±{cw} nt"
          + (", flank subtracted" if args.subtract_flank else ""))
    print(f"{'new':>4} {'file':>5}  {'sample':34}{'score':>12}{'peak height':>13}{'top 1% loci share':>19}")
    for i, pn in enumerate(top, 1):
        sc = "" if args.rank_by == "file" else f"{score[pn]:.4g}"
        print(f"{i:>4} {old_rank[pn]:>5}  {pn:34}{sc:>12}{peaks[i - 1]:>13.4g}{conc[pn]:>19.0%}")
    print(f"peak heights out of rank order: {inversions} of {pairs} pairs")
    moved = [pn for pn in top if old_rank[pn] > args.top]
    if moved:
        print(f"entered the top {args.top} relative to the file's ranking: {', '.join(moved)}")

    out = args.out or args.matrix.with_name(f"explore_{args.curve}_{args.rank_by}"
                                            + ("_flanksub" if args.subtract_flank else "") + ".pdf")
    legend = {pn: f"rank {i}, was {old_rank[pn]}" for i, pn in enumerate(top, 1)}
    ylabel = YLABELS[args.curve].format(cw=cw, sm=", smoothed" if args.sigma > 0 else "")
    if args.subtract_flank:
        ylabel += " − flank mean"
    eng.render_metaprofile(
        offsets, shown, top, legend, window, out,
        f"{args.matrix.parent.name}: {args.curve} curve, ranked by {args.rank_by}\n"
        f"top {len(top)} of {len(curves)} samples   |   n = {len(meta['loci']):,} loci",
        cw, ylabel=ylabel, hline=0.0 if args.subtract_flank else None,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
