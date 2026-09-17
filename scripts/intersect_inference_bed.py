#!/usr/bin/env python3
"""
Rank a panel of CLIP samples by how much of their binding sits at an inference BED's loci.

Workflow
  1. For each panel sample, intersect its peaks with a +/-window around every locus
     (strand-aware) and spread each peak's cDNA (score, column 5) evenly across the
     nucleotides it covers.
  2. Normalise: proportional_binding = cDNA of the sample's DISTINCT peaks inside the
     windows / the sample's peak cDNA inside --norm-bed, mitochondrial peaks excluded; each
     peak counts by the fraction of its width inside.
  3. Rank samples by central binding - mean over all loci of log1p(support within
     +/---central-window nt, per M region cDNA) - and keep the top --support-pct percent.
  4. Plot the metaprofile (first 10 of those), write sample_summary.tsv, then plot the
     heatmap and optional tSNE over the kept samples. Heatmap cells show only each locus's
     strongest central peak; the ranking and the metaprofile count every peak.

Inference loci are anchored at (start+end)//2. Panel intervals are spread across their width,
which for a 1 nt crosslink site puts the whole score at the site itself.
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
# Embed TrueType fonts in PDFs so text stays editable (e.g. in Illustrator) rather than being
# converted to Type 3 outlines.
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

CLUSTER_HUES = [
    "#393b79", # dark blue
    "#637939", # dark green
    "#e7ba52", # gold
    "#d6616b", # dark pink
    "#a55194", # purple
    "#6b6ecf", # medium blue
    "#b5cf6b", # medium green
    "#8c6d31", # brown
    "#e7969c", # light pink
    "#de9ed6", # light purple
]

DEFAULT_GENOME = "/camp/home/jonesm6/home/shared/genomes/hg38/hg38.genome"
# One fixed seed for k-means and tSNE, so a run reproduces.
RANDOM_STATE = 42
# The metaprofile draws at most this many curves, few enough that each stays distinguishable,
# while the heatmap and tSNE take every selected sample.
METAPROFILE_MAX = 10
# Mitochondrial contig names in either convention. chrM peaks never enter a sample's
# normalisation denominator: mt-rRNA is a large eCLIP background (one sample took 96% of its
# locus signal from chrM), and counting it would rank samples by how clean their library is
# rather than by where they bind.
CHRM = {"chrM", "chrMT", "MT", "M"}
PER_MILLION = 1e6
METAPROFILE_YLABEL = "Mean support per locus (per M region cDNA, smoothed)"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-x", "--xldir", required=True,
                   help="Root directory the samplesheet's file paths are resolved against.")
    p.add_argument("-b", "--bed", required=True,
                   help="Inference BED (BED6+, strand in column 6). 1 nt anchors or intervals; "
                        "intervals are anchored at their midpoint.")
    p.add_argument("-s", "--samplesheet", required=True,
                   help="TSV with columns 'file' and 'group'; 'file' is resolved relative to --xldir.")
    p.add_argument(
        "--norm-bed",
        required=True,
        help=(
            "Regions each sample's normalising cDNA is summed over (strand-aware, chrM "
            "excluded; each peak counts by the fraction of its width inside). Match it to the loci: regions_exonic.bed for an exonic locus set, "
            "regions_intronic.bed for an intronic one - both written by "
            "split_inference_bed_by_region.py."
        ),
    )
    p.add_argument("-o", "--outdir", default="results", help="Output directory (default: results/)")
    p.add_argument("--genome", default=DEFAULT_GENOME, help="Genome sizes file for bedtools slop")
    p.add_argument("--window", type=int, default=100,
                   help="Half-window size in nt around each locus (default 100)")
    p.add_argument("--gaussian-sigma", type=float, default=2.0,
                   help="Gaussian smoothing sigma for the metaprofile (default 2.0)")
    p.add_argument(
        "--support-pct",
        type=float,
        default=30.0,
        help=(
            "Keep the top P%% of panel samples by central binding for the heatmap, the "
            "tSNE and clustering (default 30). The metaprofile draws the first "
            f"{METAPROFILE_MAX} of them."
        ),
    )
    p.add_argument(
        "--central-window",
        type=int,
        default=10,
        help=(
            "Half-width (nt) of the window around nt 0 counted by the ranking score (default 10, "
            "i.e. offsets -10..+10). Every peak inside it counts toward the rank; a heatmap cell "
            "shows only the peak with the most cDNA inside it."
        ),
    )
    p.add_argument("--heatmap-scale-percentile", type=float, default=99.0,
                   help="Percentile of non-zero values mapped to the top of the colour range (default 99)")
    p.add_argument(
        "-n",
        "--n-clusters",
        type=int,
        default=None,
        help=(
            "Cluster the heatmap's loci into this many k-means groups on binarised presence. "
            "Omitted by default, when loci are ordered by total support. Passing it also "
            "writes binf_heatmap_clusters.tsv and one metaprofile per cluster, and colours "
            "the tSNE by cluster."
        ),
    )
    p.add_argument("--tsne", action="store_true",
                   help="Also write a tSNE of the heatmap's loci over the selected samples")
    p.add_argument("--tsne-perplexity", type=float, default=30.0, help="tSNE perplexity (default 30)")
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
    raise ValueError(f"Input file appears empty or has no BED rows: {path}")


def load_samplesheet_inputs(samplesheet: Path, xldir: Path) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    with open(samplesheet, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "file" not in reader.fieldnames or "group" not in reader.fieldnames:
            raise ValueError("Samplesheet must contain TSV columns: file, group")
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
    """Keep the first of each label as-is; later duplicates get _2, _3, ..."""
    counts: dict[str, int] = {}
    out: list[tuple[str, Path]] = []
    for name, path in named_paths:
        n = counts.get(name, 0) + 1
        counts[name] = n
        out.append((name if n == 1 else f"{name}_{n}", path))
    return out


def _chrom_style(path: Path, sample: int = 2000):
    """
    True if this BED is chr-prefixed, False if Ensembl-style, None if it has no data rows.

    Decided by majority over the first `sample` data rows, not the first row alone: `sort`
    is ASCII, so scaffolds (GL000009.2, KI270302.1) sort before "chr", and a chr-prefixed
    file whose first row is a scaffold looks Ensembl-style.
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
    Rewrite any panel file whose chromosome naming disagrees with the inference BED.

    The failure this prevents is silent: bedtools finds no overlap between "1" and "chr1", so
    a mismatched column reads as an RBP that binds nothing. Only mismatched files are
    rewritten, into tmpdir.
    """
    binf_chr = _chrom_style(binf_path)
    if binf_chr is None:
        return protein_sources

    fixed = []
    renamed = []
    empty = []
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
            "Chromosome naming: rewrote %d panel column(s) to match the inference BED "
            "(%s style): %s%s"
            % (len(renamed), "chr-prefixed" if binf_chr else "Ensembl",
               ", ".join(renamed[:5]), " ..." if len(renamed) > 5 else "")
        )
    if empty:
        print("WARNING: %d panel column(s) have no data rows: %s" % (len(empty), ", ".join(empty[:5])))
    return fixed


def load_binf_and_prepare_windows(binf_path: Path, window: int, genome: str, tmpdir: Path):
    """
    Anchor every inference locus at (start+end)//2 and write one +/-window BED row per
    distinct (chrom, anchor, strand).

    For a 1 nt locus the anchor is the start itself, so crosslink sites and pre-collapsed
    anchors behave exactly as before; a wider interval is measured from its middle rather
    than its start, which would otherwise smear the profile by the interval's width.

    Loci are keyed by strand as well as position. A position-only key let a + and a - locus
    at the same coordinate share one index list, so a + strand peak found by the + window
    was credited to the - locus as well.

    Returns
      binf_keys   original "chrom_start_end" per input row, in file order
      binf_index  (chrom, anchor, strand) -> row indices sharing that anchor
      windows     path to the slopped BED (7 columns, anchor last)
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
                raise ValueError(f"Inference BED has <6 columns: {line[:120]}")
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
        print(f"Inference BED: anchored {n_wide:,} interval(s) wider than 1 nt at their midpoint")
    with open(windows, "w", encoding="utf-8") as fh:
        subprocess.run(["bedtools", "slop", "-i", str(sites), "-g", genome, "-b", str(window)],
                       stdout=fh, check=True)
    return binf_keys, binf_index, windows, n_chrm


def compute_counts_for_protein(panel_bed: Path, windows_bed: Path, binf_index, window: int, n_binf: int,
                               central_window: int = 10):
    """
    counts[locus, offset] = panel peak cDNA (score, column 5) at that offset, with each peak's
    score spread evenly across its width: a peak of width w adds score/w at every nucleotide it
    covers inside the +/-window. A peak partly overlapping a window contributes only the
    overlapping part, so an 11 nt peak centred at +12 adds 4/11 of its score to offsets
    +7..+10; under a midpoint rule the central +/-10 window would have received nothing. A 1 nt
    interval (a crosslink site) keeps its whole score at one offset.

    Offsets are strand-aligned, so positive is always 5'->3' of the locus.

    Also returns locus_cdna: the cDNA of the DISTINCT peaks inside the union of all locus
    windows, each weighted by the fraction of its width inside that union. A peak inside two
    overlapping windows is added to both rows of counts but only once here, which is the right
    numerator for a proportion - 76% of the THRAP3 exonic loci have a same-strand neighbour
    within 200 nt, and summing per window inflated totals ~1.76x.

    Also returns central_max, FOR THE HEATMAP ONLY: per locus, the in-window cDNA of the peak
    with the most cDNA inside +/-central_window (score x overlap / width), so a large peak
    clipping the edge cannot beat a smaller one sitting on the locus. The ranking and the
    metaprofile use counts, where every peak is present.

    Window columns are read from the END of each intersect row, so panel files with more than
    six columns are handled.
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
    """Strand-aware merge of a region BED into BED6, so an overlap is never counted twice."""
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
    cDNA of the sample's peaks inside norm_bed (strand-aware, chrM skipped), each peak weighted
    by the fraction of its width inside the regions. norm_bed must already be strand-merged
    (merge_regions), or a peak overlapping two regions would be counted twice.
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
    Per-locus statistics over each row of counts (n_loci x L).

    Returns totals, variance, pearson_median_skew, kurtosis_excess, max_binding_offset.
    max_binding_offset is the centre of the highest-summing 5 nt sliding window, 0 for empty
    loci. On ties it takes the centre of the FIRST tied run: a lone peak fills five windows
    with the same sum, and np.argmax alone picked the leftmost, reporting it 2 nt upstream
    of where it sits.
    """
    totals = counts.sum(axis=1)
    mean = counts.mean(axis=1)
    median = np.median(counts, axis=1)
    variance = counts.var(axis=1)
    std = counts.std(axis=1)

    pearson_median_skew = np.zeros(counts.shape[0], dtype=np.float64)
    mask_std = std > 0
    pearson_median_skew[mask_std] = 3.0 * (mean[mask_std] - median[mask_std]) / std[mask_std]

    # Excess kurtosis (Fisher): mu4/sigma^4 - 3
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
    log1p, then scale to [0, 1] against a high percentile of the NON-ZERO values.

    Empty cells stay at 0, so the whole palette carries signal. The percentile is over
    non-zero values because on a sparse matrix a percentile over all cells is itself 0. The
    log1p matters: support is heavy-tailed, and scaling raw values put the median cell at 0.11.
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


def render_metaprofile(offsets, profiles, order, legend, n_loci, window, out_path, title, central_window) -> None:
    """
    One metaprofile panel. Left axis: mean normalised support per locus, every locus in the
    denominator. Right axis: the same curve times n_loci - one constant, so both axes
    describe the same pixels.
    """
    fig_h = max(9.6, 0.30 * len(order) + 2.2)
    fig, ax = plt.subplots(figsize=(13.5, fig_h))
    # 10 colours, then change linestyle when they wrap: a second non-colour channel stays
    # readable under colour-vision deficiency, which a 20-hue palette does not.
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    linestyles = ["-", "--", ":", "-."]
    for i, pn in enumerate(order):
        ax.plot(
            offsets,
            profiles[pn],
            label=f"{pn}  ({legend[pn]})",
            color=palette[i % len(palette)],
            linestyle=linestyles[(i // len(palette)) % len(linestyles)],
            linewidth=2,
        )
    ax.axvline(0, color="black", linewidth=1, alpha=0.4)
    # Mark the +/-central-window that the ranking score counts.
    for edge in (-central_window, central_window):
        ax.axvline(edge, color="red", linestyle=":", linewidth=1.2)
    ax.set_xlabel("Relative nucleotide position around inference loci (nt)")
    ax.set_ylabel(METAPROFILE_YLABEL)
    ax.set_xlim(-window, window)
    ax.set_title(title, loc="left", fontsize=10)

    ax2 = ax.twinx()
    lo, hi = ax.get_ylim()
    ax2.set_ylim(lo * n_loci, hi * n_loci)
    ax2.set_ylabel(f"Summed over {n_loci:,} loci (per M region cDNA)")

    # Figure-fraction layout: an axes-relative legend is blind to how wide ax2's tick labels
    # turn out, and wide ones pushed the y-label under the legend text.
    fig.subplots_adjust(left=0.06, right=0.50, top=1.0 - 0.5 / fig_h, bottom=0.72 / fig_h)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="center left",
               bbox_to_anchor=(0.60, 0.5), bbox_transform=fig.transFigure, borderaxespad=0.0)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_cluster_metaprofiles(protein_sources, order, scale, cluster_ids, labels, windows_bed,
                              binf_index, window, n_binf, sigma, outdir, central_window, region) -> None:
    """One normalised metaprofile per k-means cluster, over the same samples as the global one."""
    offsets = np.arange(-window, window + 1, dtype=np.int64)
    masks = {cid: labels == cid for cid in cluster_ids}
    profiles: dict[int, dict[str, np.ndarray]] = {cid: {} for cid in cluster_ids}
    shares: dict[int, dict[str, str]] = {cid: {} for cid in cluster_ids}
    wanted = set(order)
    for pn, path in protein_sources:
        if pn not in wanted:
            continue
        counts, _, _ = compute_counts_for_protein(path, windows_bed, binf_index, window, n_binf)
        for cid in cluster_ids:
            m = masks[cid]
            if not m.any():
                continue
            profiles[cid][pn] = smooth_metaprofile_gaussian(counts[m].mean(axis=0), sigma) * scale[pn]
            shares[cid][pn] = f"{float(counts[m].sum()) * scale[pn] / PER_MILLION:.2%} of region cDNA"
    for cid in cluster_ids:
        if not profiles[cid]:
            continue
        n = int(masks[cid].sum())
        out = outdir / f"metaprofile_cluster_C{cid}.pdf"
        render_metaprofile(offsets, profiles[cid], [pn for pn in order if pn in profiles[cid]],
                           shares[cid], n, window, out,
                           f"{region}: cluster C{cid} metaprofile   |   n = {n:,} loci", central_window)
        print(f"Wrote cluster metaprofile plot to: {out}")


def make_tsne(matrix, labels, cluster_to_color, out_path, perplexity) -> None:
    """tSNE of the heatmap's loci, on the same scaled matrix the heatmap draws."""
    try:
        from sklearn.manifold import TSNE
    except ImportError as e:
        raise ImportError("scikit-learn is required for --tsne.") from e
    n_distinct = int(np.unique(matrix, axis=0).shape[0])
    if n_distinct < 3:
        # Nothing to embed, and sklearn's tSNE segfaulted on a matrix whose rows were all
        # identical rather than raising.
        print(f"tSNE skipped: only {n_distinct} distinct locus profile(s) in the heatmap")
        return
    used = min(float(perplexity), float(matrix.shape[0] - 1))
    emb = TSNE(n_components=2, perplexity=used, random_state=RANDOM_STATE).fit_transform(matrix)

    plt.figure(figsize=(7, 6))
    if labels is None:
        plt.scatter(emb[:, 0], emb[:, 1], s=12, alpha=0.8, edgecolors="none")
    else:
        # Cluster id is categorical: one series per cluster, in the heatmap's colours.
        for cid in sorted(np.unique(labels)):
            m = labels == cid
            plt.scatter(emb[m, 0], emb[m, 1], s=12, alpha=0.85, color=cluster_to_color.get(int(cid)),
                        edgecolors="none", label=f"C{int(cid)}")
        plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8,
                   title="Locus cluster")
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    plt.title(f"tSNE of {matrix.shape[0]:,} loci over {matrix.shape[1]} samples")
    plt.tight_layout(rect=[0.0, 0.0, 0.85, 1.0])
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Wrote tSNE plot to: {out_path} (perplexity={used:g})")


def _fmt(x) -> str:
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.6g}"


def _mean_or_nan(values: np.ndarray) -> float:
    return float(values.mean()) if values.size else float("nan")


def main():
    args = parse_args()
    if args.n_clusters is not None and args.n_clusters < 1:
        raise ValueError("--n-clusters must be >= 1")
    if not 0 < args.support_pct <= 100:
        raise ValueError("--support-pct must be in (0, 100]")
    if args.window < 2:
        raise ValueError("--window must be >= 2 (the binding-offset window is 5 nt)")
    if not 0 <= args.central_window <= args.window:
        raise ValueError("--central-window must be between 0 and --window")

    xldir = Path(args.xldir)
    binf_path = Path(args.bed)
    norm_bed = Path(args.norm_bed)
    for label, p in [("XL directory", xldir), ("Inference BED", binf_path),
                     ("--norm-bed", norm_bed), ("Genome sizes file", Path(args.genome))]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")
    binf_style, norm_style = _chrom_style(binf_path), _chrom_style(norm_bed)
    if binf_style is not None and norm_style is not None and binf_style != norm_style:
        raise ValueError(
            f"--norm-bed uses {'chr-prefixed' if norm_style else 'Ensembl'} chromosome names but "
            f"the inference BED does not; every sample would normalise to zero."
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    protein_sources = uniquify_names(load_samplesheet_inputs(Path(args.samplesheet), xldir))
    do_cluster = args.n_clusters is not None

    with tempfile.TemporaryDirectory(prefix="intersect_binf_") as tmp:
        tmpdir = Path(tmp)
        norm_merged = merge_regions(norm_bed, tmpdir)
        binf_keys, binf_index, windows_bed, n_chrm = load_binf_and_prepare_windows(
            binf_path, args.window, args.genome, tmpdir
        )
        if n_chrm:
            print(
                f"WARNING: {n_chrm:,} inference loci are on chrM. The normalisation denominator "
                "excludes chrM, so peaks at these loci inflate proportional binding; resplit "
                "with split_inference_bed_by_region.py --drop-chrM."
            )
        protein_sources = harmonise_panel_chroms(protein_sources, binf_path, tmpdir)
        protein_names = [name for name, _ in protein_sources]
        n_binf = len(binf_keys)
        offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)
        # Offsets counted by the ranking score: 1 within +/-central-window nt of the locus, else 0.
        central = (np.abs(offsets) <= args.central_window).astype(np.float64)

        # ---- 1. per-sample counts, statistics and normalisation ----
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
            profiles[pn] = smooth_metaprofile_gaussian(counts.mean(axis=0), args.gaussian_sigma) * scale
            stats[pn] = {
                "locus_cdna": locus_cdna,
                "region_cdna": reg,
                "scale": scale,
                "prop": locus_cdna / reg if reg > 0 else float("nan"),
                "central": float(np.log1p((counts.astype(np.float64) @ central) * scale).mean())
                if reg > 0 else float("nan"),
                "total": float(totals.sum()),
                "n_sig": int(has.sum()),
                "mean_offset": _mean_or_nan(maxoff[has].astype(np.float64)),
                # Variance scales with depth, so average it on the normalised scale; skew,
                # kurtosis and offset are scale-free already.
                "mean_variance": _mean_or_nan(variance[has]) * scale**2 if reg > 0 else float("nan"),
                "mean_skew": _mean_or_nan(skew[has]),
                "mean_kurt": _mean_or_nan(kurt[has]),
            }

        no_region = [pn for pn in protein_names if stats[pn]["region_cdna"] <= 0]
        if no_region:
            print(f"WARNING: {len(no_region)} sample(s) have no peak cDNA inside --norm-bed, so their "
                  f"binding scores are NA and they rank last: {', '.join(no_region[:5])}" + (" ..." if len(no_region) > 5 else ""))

        # ---- 2. rank by central binding, keep the top --support-pct ----
        # central_binding = mean over ALL loci of log1p(support within +/-central-window nt of the
        # locus, per M region cDNA). Every peak in the window counts here; only the heatmap is
        # restricted to the strongest one. Each part answers a failure seen on the THRAP3 runs:
        #   - the central window counts binding at the locus, not binding 50-100 nt away;
        #   - dividing by region cDNA takes sequencing depth out;
        #   - log1p stops a few loci deciding the rank: raw mean support was led by K562-SSB,
        #     97% of whose support came from 22 loci;
        #   - averaging over ALL loci keeps sparse samples down (loci they miss score 0), which
        #     the plain ratio did not: K562-SUPV3L1 scored 81% on 8,594 region cDNA.
        ranked = sorted(protein_names,
                        key=lambda pn: stats[pn]["central"] if stats[pn]["region_cdna"] > 0 else -1.0,
                        reverse=True)
        rank = {pn: i + 1 for i, pn in enumerate(ranked)}
        k_sel = max(1, math.ceil(args.support_pct / 100.0 * len(ranked)))
        # Only samples with region cDNA and a central peak somewhere can be selected: either one missing gives
        # an all-zero heatmap column, and cosine distance to a zero vector is undefined. Such
        # samples already rank last, so this never changes an ordinary selection.
        selected = [pn for pn in ranked if stats[pn]["region_cdna"] > 0 and best_by[pn].any()][:k_sel]
        if not selected:
            raise ValueError("No sample has region cDNA and a central peak at these loci; nothing to plot.")
        if len(selected) < k_sel:
            print(f"Note: only {len(selected)} of the {k_sel} requested samples have region cDNA and support at these loci")
        print(f"Selected the top {len(selected)} of {len(ranked)} samples ({args.support_pct:g}%) by central binding "
              f"(±{args.central_window} nt):")
        for pn in selected:
            s = stats[pn]
            print(f"  {rank[pn]:>3}. {pn:40} central={s['central']:>8.4f}  "
                  f"mean peak support={s['total'] / n_binf:>10,.1f}  loci={s['n_sig']:>6,}")

        # ---- 3. metaprofile ----
        meta_set = selected[:METAPROFILE_MAX]
        meta_path = outdir / "metaprofile.pdf"
        render_metaprofile(
            offsets, profiles, meta_set,
            {pn: f"central binding {stats[pn]['central']:.3f}" for pn in meta_set},
            n_binf, args.window, meta_path,
            f"{binf_path.stem}: top {len(meta_set)} of {len(ranked)} samples by central binding "
            f"(±{args.central_window} nt, red dotted lines)"
            f"   |   n = {n_binf:,} loci",
            args.central_window,
        )
        print(f"Wrote metaprofile plot to: {meta_path}")

        # ---- 4. the one combined table ----
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

        # ---- 5. heatmap over the selected samples ----
        matrix = np.column_stack([best_by[pn] * stats[pn]["scale"] for pn in selected])
        row_sums = matrix.sum(axis=1)
        keep = row_sums > 0
        print(f"Heatmap: {int(keep.sum()):,} of {n_binf:,} loci have a central peak from at least one selected sample")
        if not keep.any():
            raise ValueError("No locus has support from any selected sample; nothing to plot.")
        matrix_kept = matrix[keep]
        scaled = percentile_scale(matrix_kept, args.heatmap_scale_percentile)
        row_totals = row_sums[keep]

        row_clusters = None
        cluster_ids: list[int] = []
        cluster_to_color: dict = {}
        row_colors = None
        if do_cluster:
            try:
                from sklearn.cluster import KMeans
            except ImportError as e:
                raise ImportError("scikit-learn is required for -n/--n-clusters.") from e
            binary = (matrix_kept > 0).astype(np.float64)
            # k-means cannot find more clusters than there are distinct presence patterns.
            k_fit = min(args.n_clusters, int(np.unique(binary, axis=0).shape[0]))
            km = KMeans(n_clusters=k_fit, random_state=RANDOM_STATE, n_init="auto")
            row_clusters = km.fit_predict(binary).astype(np.int64) + 1
            sort_idx = np.lexsort((-row_totals, row_clusters))
            cluster_ids = sorted(int(c) for c in np.unique(row_clusters))
            palette = CLUSTER_HUES[:len(cluster_ids)] if len(cluster_ids) <= len(CLUSTER_HUES) \
                else sns.color_palette("husl", n_colors=len(cluster_ids))
            cluster_to_color = {cid: palette[i] for i, cid in enumerate(cluster_ids)}
            row_colors = [cluster_to_color[int(c)] for c in row_clusters[sort_idx]]
            print(f"Row cluster sizes (binary k-means, k={k_fit}): "
                  f"{ {cid: int((row_clusters == cid).sum()) for cid in cluster_ids} }")
        else:
            sort_idx = np.argsort(-row_totals, kind="stable")
        display = scaled[sort_idx]

        col_linkage = linkage(pdist(scaled.T, metric="cosine"), method="average") if scaled.shape[1] > 1 else None
        cbar_label = f"log1p strongest central peak per M region cDNA\n({args.heatmap_scale_percentile:g}th pct clip)"
        n_prot = display.shape[1]
        # Transposed so sample names read horizontally on the left; loci run along the x-axis.
        heatmap_fig = sns.clustermap(
            display.T,
            row_cluster=(n_prot > 1),
            row_linkage=col_linkage,
            col_cluster=False,
            col_colors=row_colors,
            # Loci are not clustered, so the column-dendrogram axis only carries the title; the
            # default 20% of the figure height left a large gap above the heatmap.
            dendrogram_ratio=(0.2, 0.04),
            cmap="cubehelix",
            # Rows follow the sample dendrogram, not rank, so print the rank on each label.
            yticklabels=[f"{pn}  [{rank[pn]}]" for pn in selected],
            xticklabels=False,
            figsize=(11, max(5.0, 0.34 * n_prot + 2.0)),
            cbar_kws={"label": cbar_label},
            vmin=0.0,
            vmax=1.0,
            # Rasterise the cell mesh only: as vector graphics every locus x sample cell is its
            # own PDF rectangle, which is huge and slow to open. Text, dendrogram and colourbar
            # stay vector.
            rasterized=True,
        )
        if heatmap_fig.ax_col_colors is not None:
            for coll in heatmap_fig.ax_col_colors.collections:
                coll.set_rasterized(True)
        heatmap_fig.ax_heatmap.set_xlabel("Inference BED loci")
        # Move the dendrogram to the right of the heatmap so the labels can sit on the left.
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
        if do_cluster:
            _hm.legend(
                handles=[Patch(facecolor=cluster_to_color[cid], edgecolor="none", label=f"C{cid}")
                         for cid in cluster_ids],
                title="Locus cluster", loc="upper left", bbox_to_anchor=(_p_leg_x, 1.0),
                frameon=False, fontsize=8, title_fontsize=9,
            )
        heatmap_fig.ax_cbar.set_position([
            _p_hm.x0 + _p_leg_x * _p_hm.width, _p_hm.y0 + 0.05 * _p_hm.height,
            0.015, max(0.12, 0.28 * _p_hm.height),
        ])
        heatmap_fig.ax_cbar.tick_params(labelsize=7)
        heatmap_fig.ax_cbar.set_ylabel(cbar_label, fontsize=8)
        # The column-dendrogram axis is empty (loci are not clustered) and sits above any cluster
        # colour bar, so a title there never collides with the data.
        heatmap_fig.ax_col_dendrogram.set_title(
            f"{binf_path.stem}: top {len(selected)} of {len(ranked)} samples by central binding "
            f"(±{args.central_window} nt)   |   {int(keep.sum()):,} loci with a central peak",
            loc="left", fontsize=10,
        )
        heatmap_path = outdir / "binf_support_heatmap.pdf"
        heatmap_fig.savefig(heatmap_path, dpi=200, bbox_inches="tight")
        plt.close(heatmap_fig.fig)
        print(f"Wrote heatmap to: {heatmap_path}" + (" (clustered)" if do_cluster else " (loci ordered by support)"))

        # ---- 6. clustering outputs ----
        if do_cluster:
            labels_all = np.full(n_binf, -1, dtype=np.int64)
            labels_all[np.flatnonzero(keep)] = row_clusters
            cluster_tsv = outdir / "binf_heatmap_clusters.tsv"
            with open(cluster_tsv, "w", encoding="utf-8") as fout:
                fout.write("binf_chr_start_end\tchrom\tstart\tend\trow_sum_support\tpasses_heatmap_filter\theatmap_cluster\n")
                for i, key in enumerate(binf_keys):
                    chrom, start_str, end_str = key.rsplit("_", 2)
                    fout.write("\t".join([key, chrom, start_str, end_str, _fmt(float(row_sums[i])),
                                          "True" if keep[i] else "False",
                                          str(int(labels_all[i])) if keep[i] else "NA"]) + "\n")
            print(f"Wrote heatmap cluster assignments to: {cluster_tsv}")
            plot_cluster_metaprofiles(
                protein_sources, meta_set, {pn: stats[pn]["scale"] for pn in meta_set}, cluster_ids,
                labels_all, windows_bed, binf_index, args.window, n_binf, args.gaussian_sigma, outdir,
                args.central_window, binf_path.stem,
            )

        # ---- 7. tSNE ----
        if args.tsne:
            make_tsne(scaled, row_clusters, cluster_to_color, outdir / "binf_summary_tsne.png", args.tsne_perplexity)


if __name__ == "__main__":
    main()
