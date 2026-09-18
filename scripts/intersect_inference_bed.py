#!/usr/bin/env python3
"""
Rank a panel of CLIP samples by their binding at the loci of an inference BED.

Method
  1. Signal. Each locus is anchored at (start + end) // 2 and extended by --window nt on
     either side. Panel peaks are intersected with these windows on the same strand, and the
     cDNA count of each peak (BED score, column 5) is distributed uniformly across the
     nucleotides it covers. Offsets are strand-aligned, with positive values downstream.
  2. Normalisation. region_cdna is the cDNA of a sample's peaks inside --norm-bed, with each
     peak weighted by the fraction of its width inside the regions and mitochondrial peaks
     excluded. Support is expressed per million region cDNA.
  3. Metaprofile and ranking. For each sample, m(o) is the mean over all loci of
     log1p(normalised support) at offset o. central_binding is the area under m(o) within
     +/- --central-window nt. Samples are ranked by central_binding and the top --support-pct
     percent are selected. The ranking score and the plotted curve are the same statistic, so
     curve height reflects rank.
  4. Output. metaprofile.pdf (Gaussian-smoothed m(o) of the highest-ranked samples),
     sample_summary.tsv (one row per sample), binf_support_heatmap.pdf (loci x selected
     samples) and, with --tsne, binf_summary_tsne.png.
"""

import argparse
import csv
import gzip
import math
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# Embed TrueType fonts so that text in the PDFs remains editable.
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

# Seed for the tSNE embedding.
RANDOM_STATE = 42
# Maximum number of curves drawn in the metaprofile.
METAPROFILE_MAX = 10
# Mitochondrial contig names. Mitochondrial peaks are excluded from region_cdna because
# mitochondrial rRNA is a major source of background in eCLIP libraries.
CHRM = {"chrM", "chrMT", "MT", "M"}
PER_MILLION = 1e6
METAPROFILE_YLABEL = "Mean log1p support per locus\n(per M region cDNA, smoothed)"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-x", "--xldir", required=True,
                   help="Directory against which the samplesheet's file paths are resolved")
    p.add_argument("-b", "--bed", required=True,
                   help="Inference BED (BED6+, strand in column 6); each locus is anchored at its midpoint")
    p.add_argument("-s", "--samplesheet", required=True,
                   help="TSV with columns 'file' (relative to --xldir) and 'group' (sample label)")
    p.add_argument("--norm-bed", required=True,
                   help="BED6 of the regions over which each sample's normalising cDNA is summed, "
                        "e.g. exons for exonic loci (written by split_inference_bed_by_region.py)")
    p.add_argument("--genome", required=True, help="Chromosome sizes file for bedtools slop")
    p.add_argument("-o", "--outdir", default="results", help="Output directory (default: results)")
    p.add_argument("--window", type=int, default=100,
                   help="Half-width (nt) of the window around each locus (default: 100)")
    p.add_argument("--central-window", type=int, default=10,
                   help="Half-width (nt) of the central window; central_binding is the area under the "
                        "metaprofile curve within it (default: 10)")
    p.add_argument("--support-pct", type=float, default=30.0,
                   help="Percentage of top-ranked samples shown in the heatmap and tSNE (default: 30); "
                        f"the metaprofile shows the first {METAPROFILE_MAX} of them")
    p.add_argument("--gaussian-sigma", type=float, default=2.0,
                   help="Standard deviation (nt) of the Gaussian kernel used to smooth the plotted "
                        "metaprofile; the ranking uses the unsmoothed curve (default: 2)")
    p.add_argument("--heatmap-scale-percentile", type=float, default=99.0,
                   help="Percentile of non-zero heatmap values mapped to the top of the colour scale (default: 99)")
    p.add_argument("--tsne", action="store_true", help="Write a tSNE embedding of the heatmap loci")
    p.add_argument("--tsne-perplexity", type=float, default=30.0, help="tSNE perplexity (default: 30)")
    return p.parse_args()


def _open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _is_data(line: str) -> bool:
    return bool(line.strip()) and not line.startswith(("#", "track", "browser"))


def validate_bed6(path: Path) -> None:
    with _open_text_auto(path) as handle:
        for line in handle:
            if not _is_data(line):
                continue
            if len(line.rstrip("\n").split("\t")) < 6:
                raise ValueError(f"Input file is not BED6+: {path}")
            return
    raise ValueError(f"Input file has no BED records: {path}")


def load_samplesheet_inputs(samplesheet: Path, xldir: Path) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    with open(samplesheet, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "file" not in reader.fieldnames or "group" not in reader.fieldnames:
            raise ValueError("Samplesheet must contain the TSV columns: file, group")
        for row in reader:
            rel = (row.get("file") or "").strip()
            grp = (row.get("group") or "").strip()
            if not rel or not grp:
                continue
            p = (xldir / rel).resolve()
            if not p.exists():
                raise FileNotFoundError(f"Samplesheet input file not found: {p}")
            validate_bed6(p)
            inputs.append((grp, p))
    if not inputs:
        raise ValueError(f"No valid inputs found in samplesheet: {samplesheet}")
    return inputs


def uniquify_names(named_paths: list[tuple[str, Path]]) -> list[tuple[str, Path]]:
    """Make sample labels unique by appending _2, _3, ... to repeated labels."""
    counts: dict[str, int] = {}
    out: list[tuple[str, Path]] = []
    for name, path in named_paths:
        n = counts.get(name, 0) + 1
        counts[name] = n
        out.append((name if n == 1 else f"{name}_{n}", path))
    return out


def _chrom_style(path: Path, sample: int = 2000):
    """
    Return True if the BED uses chr-prefixed names, False for Ensembl-style names, and None if
    it has no records. The decision is a majority vote over the first `sample` records, because
    scaffold names sort before "chr" and a single leading record can be unrepresentative.
    """
    n_chr = n_other = 0
    with _open_text_auto(path) as fh:
        for line in fh:
            if not _is_data(line):
                continue
            if line.split("\t")[0].startswith("chr"):
                n_chr += 1
            else:
                n_other += 1
            if n_chr + n_other >= sample:
                break
    if n_chr + n_other == 0:
        return None
    return n_chr >= n_other


def harmonise_panel_chroms(protein_sources, binf_path: Path, tmpdir: Path):
    """
    Rewrite panel files whose chromosome naming differs from that of the inference BED.
    bedtools reports no overlap between "1" and "chr1", so a naming mismatch would otherwise
    yield zero signal for the affected sample without raising an error.
    """
    binf_chr = _chrom_style(binf_path)
    if binf_chr is None:
        return protein_sources

    fixed, renamed, empty = [], [], []
    for name, path in protein_sources:
        style = _chrom_style(path)
        if style is None:
            empty.append(name)
            fixed.append((name, path))
            continue
        if style == binf_chr:
            fixed.append((name, path))
            continue
        out = tmpdir / ("panel_chrfix_%s.bed" % re.sub(r"[^A-Za-z0-9_.-]", "_", name))
        with _open_text_auto(path) as fin, open(str(out), "w", encoding="utf-8") as fout:
            for line in fin:
                if not _is_data(line):
                    continue
                cols = line.rstrip("\n").split("\t")
                c = cols[0]
                if binf_chr:
                    c = c if c.startswith("chr") else "chr" + c
                    if c == "chrMT":
                        c = "chrM"
                else:
                    c = c[3:] if c.startswith("chr") else c
                cols[0] = c
                fout.write("\t".join(cols) + "\n")
        renamed.append(name)
        fixed.append((name, out))

    if renamed:
        print(
            "Chromosome naming: rewrote %d panel file(s) to match the inference BED (%s style): %s%s"
            % (len(renamed), "chr-prefixed" if binf_chr else "Ensembl",
               ", ".join(renamed[:5]), " ..." if len(renamed) > 5 else "")
        )
    if empty:
        print("WARNING: %d panel file(s) contain no records: %s" % (len(empty), ", ".join(empty[:5])))
    return fixed


def load_binf_and_prepare_windows(binf_path: Path, window: int, genome: str, tmpdir: Path):
    """
    Anchor each inference locus at (start + end) // 2 and write one +/-window interval per
    distinct (chromosome, anchor, strand).

    Returns
      binf_keys   "chrom_start_end" of each input record, in file order
      binf_index  (chromosome, anchor, strand) -> indices of the records sharing that anchor
      windows     path of the window BED (seven columns, anchor in the last)
      n_chrm      number of loci on a mitochondrial contig
    """
    sites = tmpdir / "binf_site.bed"
    windows = tmpdir / "binf_windows.bed"
    binf_keys: list[str] = []
    binf_index: dict[tuple[str, int, str], list[int]] = {}
    n_wide = n_chrm = 0

    with _open_text_auto(binf_path) as fin, open(sites, "w", encoding="utf-8") as fout:
        for line in fin:
            if not _is_data(line):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                raise ValueError(f"Inference BED record has fewer than 6 columns: {line[:120]}")
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            anchor = (start + end) // 2
            n_wide += (end - start) > 1
            n_chrm += chrom in CHRM
            key = (chrom, anchor, strand)
            if key not in binf_index:
                binf_index[key] = []
                fout.write("\t".join([chrom, str(anchor), str(anchor + 1), ".", ".", strand, str(anchor)]) + "\n")
            binf_index[key].append(len(binf_keys))
            binf_keys.append(f"{chrom}_{start}_{end}")

    if n_wide:
        print(f"Inference BED: {n_wide:,} interval(s) wider than 1 nt were anchored at their midpoint")
    with open(windows, "w", encoding="utf-8") as fh:
        subprocess.run(["bedtools", "slop", "-i", str(sites), "-g", genome, "-b", str(window)],
                       stdout=fh, check=True)
    return binf_keys, binf_index, windows, n_chrm


def compute_counts_for_protein(panel_bed: Path, windows_bed: Path, binf_index, window: int, n_binf: int,
                               central_window: int = 10):
    """
    Build the locus x offset support matrix of one panel sample.

    counts[locus, offset] is the peak cDNA at that offset. The score of a peak of width w is
    distributed uniformly, score / w per nucleotide, so a peak that partly overlaps a window
    contributes only the overlapping fraction, and a 1 nt interval retains its full score at
    its own position. Offsets are strand-aligned.

    Returns
      counts       float32 array, n_loci x (2 * window + 1)
      locus_cdna   cDNA of the distinct peaks inside the union of all locus windows, each
                   weighted by the fraction of its width inside that union; a peak that lies
                   in two overlapping windows is counted once
      central_max  per locus, the largest single-peak contribution within +/-central_window
                   (score x overlap / width); used for the heatmap

    Window columns are read from the end of each intersect record, so panel files with more
    than six columns are supported.
    """
    counts = np.zeros((n_binf, 2 * window + 1), dtype=np.float32)
    central_max = np.zeros(n_binf, dtype=np.float64)
    covered: dict[tuple[str, int, int, str], list] = {}
    cmd = ["bedtools", "intersect", "-a", str(panel_bed), "-b", str(windows_bed), "-s", "-wa", "-wb"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    assert proc.stdout is not None
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        w_chrom, w_strand, anchor = f[-7], f[-2], int(f[-1])
        idx_list = binf_index.get((w_chrom, anchor, w_strand))
        if not idx_list:
            continue
        start, end = int(f[1]), int(f[2])
        if end <= start:
            end = start + 1
        width = end - start
        lo, hi = max(start, anchor - window), min(end - 1, anchor + window)
        if lo > hi:
            continue
        a, b = (anchor - hi, anchor - lo) if w_strand == "-" else (lo - anchor, hi - anchor)
        score = float(f[4])
        counts[idx_list, a + window:b + window + 1] += score / width
        clo, chi = max(start, anchor - central_window), min(end - 1, anchor + central_window)
        if clo <= chi:
            central_max[idx_list] = np.maximum(central_max[idx_list], score * (chi - clo + 1) / width)
        entry = covered.get((f[0], start, end, f[5]))
        if entry is None:
            covered[(f[0], start, end, f[5])] = [score, width, [(lo, hi)]]
        else:
            entry[2].append((lo, hi))
    _, stderr = proc.communicate()
    if proc.returncode not in (0, None):
        raise RuntimeError(f"bedtools intersect failed (code={proc.returncode}): {stderr[:500]}")

    locus_cdna = 0.0
    for score, width, spans in covered.values():
        spans.sort()
        inside = 0
        cur_lo, cur_hi = spans[0]
        for lo, hi in spans[1:]:
            if lo > cur_hi + 1:
                inside += cur_hi - cur_lo + 1
                cur_lo, cur_hi = lo, hi
            else:
                cur_hi = max(cur_hi, hi)
        inside += cur_hi - cur_lo + 1
        locus_cdna += score * inside / width
    return counts, locus_cdna, central_max


def merge_regions(bed: Path, tmpdir: Path) -> Path:
    """Merge a region BED on each strand into BED6, so that overlapping regions are counted once."""
    rows = []
    with _open_text_auto(bed) as fin:
        for line in fin:
            if not _is_data(line):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                raise ValueError(f"--norm-bed must be BED6 with strand in column 6: {line[:120]}")
            rows.append((c[0], int(c[1]), int(c[2]), c[5]))
    rows.sort(key=lambda r: (r[0], r[1]))
    srt = tmpdir / "norm_sorted.bed"
    with open(srt, "w", encoding="utf-8") as fout:
        fout.writelines(f"{c}\t{a}\t{b}\t.\t.\t{st}\n" for c, a, b, st in rows)
    r = subprocess.run(["bedtools", "merge", "-s", "-c", "6", "-o", "distinct", "-i", str(srt)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools merge (normalisation regions) failed: {r.stderr[:500]}")
    merged = tmpdir / "norm_merged.bed"
    with open(merged, "w", encoding="utf-8") as fout:
        for ln in r.stdout.splitlines():
            if ln:
                c = ln.split("\t")
                fout.write(f"{c[0]}\t{c[1]}\t{c[2]}\t.\t.\t{c[3]}\n")
    return merged


def region_cdna(panel_bed: Path, norm_bed: Path, tmpdir: Path) -> float:
    """
    cDNA of a sample's peaks inside norm_bed on the same strand, excluding mitochondrial peaks.
    Each peak is weighted by the fraction of its width inside the regions. norm_bed must be
    strand-merged (merge_regions) so that no overlap is counted twice.
    """
    peaks = tmpdir / "norm_peaks.bed"
    info: list[tuple[float, int]] = []
    with _open_text_auto(panel_bed) as fin, open(peaks, "w", encoding="utf-8") as fout:
        for line in fin:
            if not _is_data(line):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] in CHRM:
                continue
            start, end = int(c[1]), int(c[2])
            if end <= start:
                end = start + 1
            fout.write(f"{c[0]}\t{start}\t{end}\t{len(info)}\t{c[4]}\t{c[5]}\n")
            info.append((float(c[4]), end - start))
    r = subprocess.run(["bedtools", "intersect", "-s", "-wo", "-a", str(peaks), "-b", str(norm_bed)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools intersect (normalisation) failed: {r.stderr[:500]}")
    overlap: dict[int, int] = {}
    for ln in r.stdout.splitlines():
        if ln:
            f = ln.split("\t")
            i = int(f[3])
            overlap[i] = overlap.get(i, 0) + int(f[-1])
    return float(sum(info[i][0] * min(n, info[i][1]) / info[i][1] for i, n in overlap.items()))


def compute_summary_stats(counts: np.ndarray, window: int):
    """
    Per-locus statistics of the support vector (rows of counts).

    Returns totals, variance, pearson_median_skew, kurtosis_excess and max_binding_offset.
    max_binding_offset is the centre of the 5 nt sliding window with the highest sum, and 0
    for loci without signal. When several adjacent windows tie, the centre of the first tied
    run is reported, so that an isolated peak is assigned to its own position.
    """
    totals = counts.sum(axis=1)
    mean = counts.mean(axis=1)
    median = np.median(counts, axis=1)
    variance = counts.var(axis=1)
    std = counts.std(axis=1)

    pearson_median_skew = np.zeros(counts.shape[0], dtype=np.float64)
    mask_std = std > 0
    pearson_median_skew[mask_std] = 3.0 * (mean[mask_std] - median[mask_std]) / std[mask_std]

    # Excess kurtosis (Fisher): mu4 / sigma^4 - 3
    centered = counts - mean[:, None]
    mu2 = np.mean(centered**2, axis=1)
    mu4 = np.mean(centered**4, axis=1)
    kurtosis_excess = np.zeros(counts.shape[0], dtype=np.float64)
    mask_var = mu2 > 0
    kurtosis_excess[mask_var] = mu4[mask_var] / (mu2[mask_var] ** 2) - 3.0

    if counts.shape[1] < 5:
        max_binding_offset = np.zeros(counts.shape[0], dtype=np.int64)
    else:
        ws = np.lib.stride_tricks.sliding_window_view(counts, 5, axis=1).sum(axis=2)
        is_max = ws == ws.max(axis=1, keepdims=True)
        first = is_max.argmax(axis=1)
        col = np.arange(ws.shape[1])
        breaks = (~is_max) & (col[None, :] > first[:, None])
        run_end = np.where(breaks.any(axis=1), breaks.argmax(axis=1), ws.shape[1])
        max_binding_offset = ((first + run_end - 1) // 2 + 2 - window).astype(np.int64)
    max_binding_offset[totals == 0] = 0

    return totals, variance, pearson_median_skew, kurtosis_excess, max_binding_offset


def smooth_metaprofile_gaussian(meta_counts: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    radius = max(1, int(round(4.0 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()
    return np.convolve(meta_counts, kernel, mode="same")


def percentile_scale(matrix: np.ndarray, pct: float) -> np.ndarray:
    """
    Apply log1p and scale to [0, 1] against the pct-th percentile of the non-zero values, with
    values above that percentile clipped to 1. Empty cells remain at 0. The percentile is taken
    over non-zero values because the matrix is sparse.
    """
    dense = np.log1p(np.clip(matrix.astype(np.float64), 0.0, None))
    nonzero = dense[dense > 0]
    if nonzero.size == 0:
        return np.zeros_like(dense)
    hi = float(np.percentile(nonzero, pct))
    if hi <= 0:
        hi = float(nonzero.max())
    if hi <= 0:
        return np.zeros_like(dense)
    return np.clip(dense / hi, 0.0, 1.0)


def log_mean_curve(counts, scale):
    """
    Metaprofile curve m(o): at each offset, the mean over all loci of log1p(support x scale),
    where scale = 1e6 / region_cdna. Unsmoothed.

    The log1p transform is applied per locus before averaging so that the curve reflects the
    breadth of binding across loci rather than being dominated by a small number of strongly
    bound loci. Loci without signal contribute 0, and the normalisation removes the dependence
    on sequencing depth. central_binding is the area under this curve within the central window.
    """
    return np.log1p(counts.astype(np.float64) * scale).mean(axis=0)


def render_metaprofile(offsets, profiles, order, legend, window, out_path, title, central_window) -> None:
    """Draw the metaprofile: square plot area, 12 pt text, legend to the right."""
    with plt.rc_context({"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12,
                         "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12}):
        fig = plt.figure(figsize=(15.0, 8.6))
        # Plot area of 7 x 7 in, given as [left, bottom, width, height] in figure fractions.
        ax = fig.add_axes([1.2 / 15.0, 0.9 / 8.6, 7.0 / 15.0, 7.0 / 8.6])
        ax.set_box_aspect(1)
        # The colour cycle has 10 entries; the line style changes each time it repeats.
        palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        linestyles = ["-", "--", ":", "-."]
        for i, pn in enumerate(order):
            ax.plot(offsets, profiles[pn], label=f"{pn}  ({legend[pn]})",
                    color=palette[i % len(palette)],
                    linestyle=linestyles[(i // len(palette)) % len(linestyles)], linewidth=2)
        ax.axvline(0, color="black", linewidth=1, alpha=0.4)
        # Limits of the central window used for ranking.
        for edge in (-central_window, central_window):
            ax.axvline(edge, color="red", linestyle=":", linewidth=1.2)
        ax.set_xlabel("Relative nucleotide position around inference loci (nt)")
        ax.set_ylabel(METAPROFILE_YLABEL)
        ax.set_xlim(-window, window)
        ax.set_title(title, loc="left")
        fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="center left",
                   bbox_to_anchor=(8.6 / 15.0, 0.5), bbox_transform=fig.transFigure, borderaxespad=0.0)
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def make_tsne(matrix, out_path, perplexity) -> None:
    """tSNE embedding of the heatmap loci, computed on the scaled heatmap matrix."""
    try:
        from sklearn.manifold import TSNE
    except ImportError as e:
        raise ImportError("scikit-learn is required for --tsne.") from e
    n_distinct = int(np.unique(matrix, axis=0).shape[0])
    if n_distinct < 3:
        print(f"tSNE skipped: the heatmap contains only {n_distinct} distinct locus profile(s)")
        return
    used = min(float(perplexity), float(matrix.shape[0] - 1))
    emb = TSNE(n_components=2, perplexity=used, random_state=RANDOM_STATE).fit_transform(matrix)

    plt.figure(figsize=(7, 6))
    plt.scatter(emb[:, 0], emb[:, 1], s=12, alpha=0.8, edgecolors="none")
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    plt.title(f"tSNE of {matrix.shape[0]:,} loci over {matrix.shape[1]} samples")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Wrote tSNE plot to: {out_path} (perplexity={used:g})")


def _fmt(x) -> str:
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.6g}"


def _mean_or_nan(values: np.ndarray) -> float:
    return float(values.mean()) if values.size else float("nan")


def main():
    args = parse_args()
    if not 0 < args.support_pct <= 100:
        raise ValueError("--support-pct must be in (0, 100]")
    if args.window < 2:
        raise ValueError("--window must be >= 2")
    if not 0 <= args.central_window <= args.window:
        raise ValueError("--central-window must be between 0 and --window")

    xldir = Path(args.xldir)
    binf_path = Path(args.bed)
    norm_bed = Path(args.norm_bed)
    for label, p in [("--xldir", xldir), ("Inference BED", binf_path),
                     ("--norm-bed", norm_bed), ("Genome sizes file", Path(args.genome))]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")
    binf_style, norm_style = _chrom_style(binf_path), _chrom_style(norm_bed)
    if binf_style is not None and norm_style is not None and binf_style != norm_style:
        raise ValueError(
            f"--norm-bed uses {'chr-prefixed' if norm_style else 'Ensembl'} chromosome names but "
            "the inference BED does not; region_cdna would be zero for every sample."
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    protein_sources = uniquify_names(load_samplesheet_inputs(Path(args.samplesheet), xldir))

    with tempfile.TemporaryDirectory(prefix="intersect_binf_") as tmp:
        tmpdir = Path(tmp)
        norm_merged = merge_regions(norm_bed, tmpdir)
        binf_keys, binf_index, windows_bed, n_chrm = load_binf_and_prepare_windows(
            binf_path, args.window, args.genome, tmpdir
        )
        if n_chrm:
            print(
                f"WARNING: {n_chrm:,} inference loci are on a mitochondrial contig. region_cdna excludes "
                "mitochondrial peaks, so binding at these loci is over-weighted; remove them with "
                "split_inference_bed_by_region.py --drop-chrM."
            )
        protein_sources = harmonise_panel_chroms(protein_sources, binf_path, tmpdir)
        protein_names = [name for name, _ in protein_sources]
        n_binf = len(binf_keys)
        offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)
        in_central = np.abs(offsets) <= args.central_window

        # 1. Support matrix, normalisation and per-sample statistics.
        best_by: dict[str, np.ndarray] = {}
        profiles: dict[str, np.ndarray] = {}
        stats: dict[str, dict] = {}
        for pn, path in protein_sources:
            counts, locus_cdna, central_max = compute_counts_for_protein(
                path, windows_bed, binf_index, args.window, n_binf, args.central_window)
            totals, variance, skew, kurt, maxoff = compute_summary_stats(counts, args.window)
            reg = region_cdna(path, norm_merged, tmpdir)
            scale = PER_MILLION / reg if reg > 0 else 0.0
            has = totals > 0
            best_by[pn] = central_max
            curve = log_mean_curve(counts, scale)
            profiles[pn] = smooth_metaprofile_gaussian(curve, args.gaussian_sigma)
            stats[pn] = {
                "locus_cdna": locus_cdna,
                "region_cdna": reg,
                "scale": scale,
                "prop": locus_cdna / reg if reg > 0 else float("nan"),
                "central": float(curve[in_central].sum()) if reg > 0 else float("nan"),
                "total": float(totals.sum()),
                "n_sig": int(has.sum()),
                "mean_offset": _mean_or_nan(maxoff[has].astype(np.float64)),
                # The variance depends on sequencing depth and is therefore reported on the
                # normalised scale; skew, kurtosis and offset are scale-free.
                "mean_variance": _mean_or_nan(variance[has]) * scale**2 if reg > 0 else float("nan"),
                "mean_skew": _mean_or_nan(skew[has]),
                "mean_kurt": _mean_or_nan(kurt[has]),
            }

        no_region = [pn for pn in protein_names if stats[pn]["region_cdna"] <= 0]
        if no_region:
            print(f"WARNING: {len(no_region)} sample(s) have no peak cDNA inside --norm-bed; their scores "
                  f"are NA and they are ranked last: {', '.join(no_region[:5])}"
                  + (" ..." if len(no_region) > 5 else ""))

        # 2. Rank by central_binding and select the top --support-pct percent.
        ranked = sorted(protein_names,
                        key=lambda pn: stats[pn]["central"] if stats[pn]["region_cdna"] > 0 else -1.0,
                        reverse=True)
        rank = {pn: i + 1 for i, pn in enumerate(ranked)}
        k_sel = max(1, math.ceil(args.support_pct / 100.0 * len(ranked)))
        # A sample without region cDNA, or without any peak in a central window, would give an
        # all-zero heatmap column, for which the cosine distance is undefined.
        selected = [pn for pn in ranked if stats[pn]["region_cdna"] > 0 and best_by[pn].any()][:k_sel]
        if not selected:
            raise ValueError("No sample has region cDNA and a peak in a central window; nothing to plot.")
        if len(selected) < k_sel:
            print(f"Note: only {len(selected)} of the {k_sel} requested samples have region cDNA and a "
                  "peak in a central window")
        print(f"Selected the top {len(selected)} of {len(ranked)} samples ({args.support_pct:g}%) by central "
              f"binding (±{args.central_window} nt):")
        for pn in selected:
            s = stats[pn]
            print(f"  {rank[pn]:>3}. {pn:40} central binding={s['central']:>8.4f}  "
                  f"mean peak support={s['total'] / n_binf:>10,.1f}  loci={s['n_sig']:>6,}")

        # 3. Metaprofile.
        meta_set = selected[:METAPROFILE_MAX]
        legend = {pn: f"rank {rank[pn]}, central binding {stats[pn]['central']:.3f}" for pn in meta_set}
        meta_path = outdir / "metaprofile.pdf"
        render_metaprofile(
            offsets, profiles, meta_set, legend, args.window, meta_path,
            f"{binf_path.stem}\ntop {len(meta_set)} of {len(ranked)} samples by central binding "
            f"(±{args.central_window} nt, red dotted lines)   |   n = {n_binf:,} loci",
            args.central_window,
        )
        print(f"Wrote metaprofile plot to: {meta_path}")

        # 4. Sample summary table.
        table_path = outdir / "sample_summary.tsv"
        cols = ["sample", "rank", "selected", "central_binding", "proportional_binding", "locus_cdna", "region_cdna",
                "mean_peak_support", "total_peak_support", "loci_with_signal", "frac_loci_with_signal",
                "mean_binding_offset", "mean_variance", "mean_pearson_skew", "mean_kurtosis"]
        chosen = set(selected)
        with open(table_path, "w", encoding="utf-8") as fout:
            fout.write("\t".join(cols) + "\n")
            for pn in ranked:
                s = stats[pn]
                fout.write("\t".join([
                    pn, str(rank[pn]), "True" if pn in chosen else "False", _fmt(s["central"]),
                    _fmt(s["prop"]), _fmt(s["locus_cdna"]), _fmt(s["region_cdna"]),
                    _fmt(s["total"] / n_binf), _fmt(s["total"]), str(s["n_sig"]), _fmt(s["n_sig"] / n_binf),
                    _fmt(s["mean_offset"]), _fmt(s["mean_variance"]), _fmt(s["mean_skew"]), _fmt(s["mean_kurt"]),
                ]) + "\n")
        print(f"Wrote sample summary to: {table_path}")

        # 5. Heatmap of the selected samples. Each cell is the strongest single peak within the
        #    central window of that locus, per million region cDNA.
        matrix = np.column_stack([best_by[pn] * stats[pn]["scale"] for pn in selected])
        row_sums = matrix.sum(axis=1)
        keep = row_sums > 0
        print(f"Heatmap: {int(keep.sum()):,} of {n_binf:,} loci have a central peak from at least one selected sample")
        if not keep.any():
            raise ValueError("No locus has a central peak from any selected sample; nothing to plot.")
        scaled = percentile_scale(matrix[keep], args.heatmap_scale_percentile)
        display = scaled[np.argsort(-row_sums[keep], kind="stable")]

        col_linkage = linkage(pdist(scaled.T, metric="cosine"), method="average") if scaled.shape[1] > 1 else None
        cbar_label = f"log1p strongest central peak per M region cDNA\n({args.heatmap_scale_percentile:g}th pct clip)"
        n_prot = display.shape[1]
        # The matrix is transposed so that sample names are read horizontally; loci run along
        # the x-axis. Only the cell mesh is rasterised; text, dendrogram and colour bar remain
        # vector graphics.
        heatmap_fig = sns.clustermap(
            display.T,
            row_cluster=(n_prot > 1),
            row_linkage=col_linkage,
            col_cluster=False,
            # Loci are not clustered, so the column-dendrogram axis only carries the title.
            dendrogram_ratio=(0.2, 0.04),
            cmap="cubehelix",
            # Rows follow the sample dendrogram, so each label carries the rank of the sample.
            yticklabels=[f"{pn}  [{rank[pn]}]" for pn in selected],
            xticklabels=False,
            figsize=(11, max(5.0, 0.34 * n_prot + 2.0)),
            cbar_kws={"label": cbar_label},
            vmin=0.0,
            vmax=1.0,
            rasterized=True,
        )
        heatmap_fig.ax_heatmap.set_xlabel("Inference BED loci")
        # Place the dendrogram to the right of the heatmap and the sample labels to the left.
        _hm = heatmap_fig.ax_heatmap
        _rd = heatmap_fig.ax_row_dendrogram
        _p_hm, _p_rd = _hm.get_position(), _rd.get_position()
        _rd.set_position([_p_hm.x1 + 0.012, _p_rd.y0, _p_rd.width, _p_rd.height])
        _rd.invert_xaxis()
        _hm.yaxis.tick_left()
        _hm.yaxis.set_label_position("left")
        _hm.set_ylabel("Samples  [n] = rank by central binding")
        plt.setp(_hm.get_yticklabels(), rotation=0, fontsize=8)
        _p_leg_x = 1.0 + (_p_rd.width / max(_p_hm.width, 1e-9)) + 0.05
        heatmap_fig.ax_cbar.set_position([
            _p_hm.x0 + _p_leg_x * _p_hm.width, _p_hm.y0 + 0.05 * _p_hm.height,
            0.015, max(0.12, 0.28 * _p_hm.height),
        ])
        heatmap_fig.ax_cbar.tick_params(labelsize=7)
        heatmap_fig.ax_cbar.set_ylabel(cbar_label, fontsize=8)
        heatmap_fig.ax_col_dendrogram.set_title(
            f"{binf_path.stem}: top {len(selected)} of {len(ranked)} samples by central binding "
            f"(±{args.central_window} nt)   |   {int(keep.sum()):,} loci with a central peak",
            loc="left", fontsize=10,
        )
        heatmap_path = outdir / "binf_support_heatmap.pdf"
        heatmap_fig.savefig(heatmap_path, dpi=200, bbox_inches="tight")
        plt.close(heatmap_fig.fig)
        print(f"Wrote heatmap to: {heatmap_path}")

        # 6. tSNE of the heatmap loci.
        if args.tsne:
            make_tsne(scaled, outdir / "binf_summary_tsne.png", args.tsne_perplexity)


if __name__ == "__main__":
    main()
