import argparse
import csv
import gzip
import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

DEFAULT_GENOME = "/camp/home/jonesm6/home/shared/genomes/hg38/hg38.genome"
HEATMAP_TOP_PROTEINS = 100
METAPROFILE_TOP_PROTEINS = 15

def parse_args():
    p = argparse.ArgumentParser(description="Intersect inference BED with multiple crosslink types and summarize")
    p.add_argument(
        "-x",
        "--xldir",
        required=True,
        help="Base XL directory. Used for protein subdirectories (merge mode) or direct BED inputs (skip-merge mode).",
    )
    p.add_argument(
        "-b",
        "--bed",
        required=True,
        help="Inference BED (must have at least 6 columns; strand in column 6).",
    )
    p.add_argument("--window", type=int, default=100, help="Half-window size in bp around each inference locus")
    p.add_argument("--gaussian-sigma", type=float, default=2.0, help="Gaussian smoothing sigma for metaprofile")
    p.add_argument(
        "--cluster-metaprofiles",
        action="store_true",
        help="If set, write one metaprofile plot per heatmap row cluster.",
    )
    p.add_argument(
        "--n-clusters",
        type=int,
        default=4,
        help="Number of row clusters to extract from hierarchical heatmap clustering (default: 4).",
    )
    p.add_argument(
        "--cluster-top-proteins",
        type=int,
        default=HEATMAP_TOP_PROTEINS,
        help=(
            "Use top K XL groups by summed total_overlaps across all loci for heatmap/clustering/tSNE "
            f"(default: {HEATMAP_TOP_PROTEINS})."
        ),
    )
    p.add_argument(
        "--metaprofile-top-proteins",
        type=int,
        default=METAPROFILE_TOP_PROTEINS,
        help=(
            "Use top K proteins for global and per-cluster metaprofile plotting "
            f"(default: {METAPROFILE_TOP_PROTEINS})."
        ),
    )
    p.add_argument(
        "-i",
        "--inspect-protein",
        default=None,
        help="Optional protein name to generate per-nucleotide binf heatmap (rows=binf, cols=-window..+window).",
    )
    p.add_argument("--genome", default=DEFAULT_GENOME, help="Genome sizes file for bedtools slop")
    p.add_argument(
        "--table",
        action="store_true",
        help="If set, write per-locus summary TSV.",
    )
    p.add_argument(
        "--tsne",
        action="store_true",
        help="If set, generate a tSNE plot from <outdir>/binf_summary.tsv total-overlap columns.",
    )
    p.add_argument(
        "--tsne-perplexity",
        type=float,
        default=30.0,
        help="tSNE perplexity (will be clipped to valid range for sample count).",
    )
    p.add_argument(
        "--tsne-random-state",
        type=int,
        default=42,
        help="Random seed for tSNE reproducibility.",
    )
    p.add_argument(
        "-o",
        "--outdir",
        default="results",
        help="Output directory for metaprofile and (optional) TSV (default: results/ in working directory)",
    )
    p.add_argument(
        "--skip-merge",
        action="store_true",
        help="Skip per-protein merge and use direct BED/BED.GZ inputs.",
    )
    p.add_argument(
        "-s",
        "--samplesheet",
        default=None,
        help=(
            "Optional TSV with columns 'file' and 'group'. "
            "'file' is resolved relative to --xldir and 'group' is used as the protein label."
        ),
    )
    return p.parse_args()


def list_protein_dirs(xldir: Path) -> list[Path]:
    """
    Returns per-protein directories under xldir.
    If xldir itself contains *genome.xl.bed and no protein subdirectories match, treat xldir as a single protein.
    """
    protein_dirs: list[Path] = []
    for child in sorted(xldir.iterdir()):
        if not child.is_dir():
            continue
        if list(child.glob("*genome.xl.bed")):
            protein_dirs.append(child)

    if protein_dirs:
        return protein_dirs

    # Fallback for the legacy layout where xldir is already the protein directory.
    if list(xldir.glob("*genome.xl.bed")):
        return [xldir]

    raise FileNotFoundError(f"No protein crosslink beds (*genome.xl.bed) found under {xldir}")


def _open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def validate_bed6(path: Path) -> None:
    with _open_text_auto(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
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


def discover_bed_inputs_from_xldir(xldir: Path) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    for p in sorted(xldir.iterdir()):
        if not p.is_file():
            continue
        name = p.name.lower()
        if not (name.endswith(".bed") or name.endswith(".bed.gz")):
            continue
        validate_bed6(p)
        # Derive group from basename without .bed/.bed.gz
        stem = p.name
        if stem.endswith(".bed.gz"):
            stem = stem[: -len(".bed.gz")]
        elif stem.endswith(".bed"):
            stem = stem[: -len(".bed")]
        inputs.append((stem, p))
    if not inputs:
        raise FileNotFoundError(f"No BED/BED.GZ files found directly in {xldir}")
    return inputs


def uniquify_names(named_paths: list[tuple[str, Path]]) -> list[tuple[str, Path]]:
    """
    Ensure unique protein labels across multiple xldirs.
    Keeps first label as-is; subsequent duplicates get _2, _3, ...
    """
    counts: dict[str, int] = {}
    out: list[tuple[str, Path]] = []
    for name, path in named_paths:
        n = counts.get(name, 0) + 1
        counts[name] = n
        if n == 1:
            out.append((name, path))
        else:
            out.append((f"{name}_{n}", path))
    return out


def merge_protein_crosslinks(protein_dir: Path, merged_dir: Path) -> Path:
    """
    Merge all *genome.xl.bed files in a protein directory into a single BED6:
      chr start end . fileSupport strand
    fileSupport = number of distinct XL files where the exact site+strand occurs.
    Chromosomes are normalized to chr-prefixed naming.
    """
    genome_xl = sorted(protein_dir.glob("*genome.xl.bed"))
    if not genome_xl:
        raise FileNotFoundError(f"No *genome.xl.bed found in {protein_dir}")

    merged_dir.mkdir(parents=True, exist_ok=True)
    protein_name = protein_dir.name
    out_path = merged_dir / f"{protein_name}_merged.bed"
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"Using existing merged BED for {protein_name}: {out_path}")
        return out_path

    # Tag each XL record with its source filename, then group by genomic site+strand
    # and count distinct source files (binary support across files).
    merge_cmd = (
        "for f in *.genome.xl.bed; do "
        "awk -v OFS='\\t' -v sid=\"$f\" '{print $1,$2,$3,$4,$5,$6,sid}' \"$f\"; "
        "done | "
        "sort -k1,1 -k2,2n -k3,3n -k6,6 -k7,7 | "
        "bedtools groupby -i stdin -g 1,2,3,6 -c 7 -o count_distinct | "
        "awk 'BEGIN{OFS=\"\\t\"} "
        "{chrom=$1; if (chrom !~ /^chr/) chrom=\"chr\" chrom; print chrom, $2, $3, \".\", $5, $4}' "
        f"> \"{out_path}\""
    )
    subprocess.run(merge_cmd, cwd=str(protein_dir), shell=True, check=True)
    return out_path


def load_binf_and_prepare_windows(binf_path: Path, window: int, genome: str, tmpdir: Path):
    """
    Writes:
      - binf_site.bed: inference BED with extra columns containing site start/end
      - binf_windows.bed: bedtools slop output of binf_site.bed (slopped intervals kept with extra columns)

    Returns:
      binf_keys: list of chr_start_end strings in stable file order
      binf_index: dict mapping (chr, start, end) -> list of row indices (duplicates allowed)
      binf_windows_path: path to slopped inference BED
    """
    binf_site = tmpdir / "binf_site.bed"
    binf_windows = tmpdir / "binf_windows.bed"

    binf_keys: list[str] = []
    # Allows duplicate loci (same chr/start/end can appear multiple times).
    # Maps (chr, start, end) -> list of row indices that share these coordinates.
    binf_index: dict[tuple[str, int, int], list[int]] = {}
    binf_idx = 0

    with open(binf_path, "r", encoding="utf-8") as fin, open(binf_site, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                raise ValueError(f"Inference BED has <6 columns: {line[:120]}")
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]
            # columns 4-5 are preserved (name and spliceAI score), plus we add site start/end
            key = f"{chrom}_{start}_{end}"
            binf_keys.append(key)
            binf_index.setdefault((chrom, start, end), []).append(binf_idx)
            # Ensure strand stays in column 6.
            fout.write("\t".join([chrom, str(start), str(end), cols[3], cols[4], strand, str(start), str(end)]) + "\n")
            binf_idx += 1

    subprocess.run(
        ["bedtools", "slop", "-i", str(binf_site), "-g", genome, "-b", str(window)],
        stdout=open(binf_windows, "w", encoding="utf-8"),
        check=True,
    )

    return binf_keys, binf_index, binf_windows


def compute_counts_for_protein(
    merged_xl_bed: Path,
    binf_windows_bed: Path,
    binf_index: dict[tuple[str, int, int], list[int]],
    window: int,
    n_binf: int,
):
    """
    Build a matrix counts[binf_idx, offset_bin] where:
      offset_bin corresponds to offset in [-window..+window]
      value = summed crosslink score (merged_xl column 5) at xl_start for that offset
    Uses strand-aware overlaps (-s) and assigns offsets with strand flipped so that offsets always match the inference BED's 5'->3' direction.
    """
    L = 2 * window + 1
    counts = np.zeros((n_binf, L), dtype=np.float32)

    # intersect output columns:
    # merged_xl: 6 columns
    # binf_windows: 8 columns (chr,start,end,name,score,strand,site_start,site_end)
    # plus overlap length (-wo)
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(merged_xl_bed),
        "-b",
        str(binf_windows_bed),
        "-s",
        "-wo",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)

    assert proc.stdout is not None
    for line in proc.stdout:
        if not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 15:
            continue

        xl_start = int(cols[1])
        xl_score = float(cols[4])
        binf_chr = cols[6]
        binf_strand = cols[11]
        binf_site_start = int(cols[12])
        binf_site_end = int(cols[13])
        idx_list = binf_index.get((binf_chr, binf_site_start, binf_site_end))
        if not idx_list:
            # Chromosome naming mismatch should not happen, but skip if it does.
            continue

        offset = xl_start - binf_site_start
        if binf_strand == "-":
            offset = -offset

        if -window <= offset <= window:
            bin_idx = offset + window
            for idx in idx_list:
                counts[idx, bin_idx] += xl_score

    _, stderr = proc.communicate()
    if proc.returncode not in (0, None):
        raise RuntimeError(f"bedtools intersect failed (code={proc.returncode}): {stderr[:500]}")

    return counts


def compute_summary_stats(counts: np.ndarray, window: int):
    """
    counts shape: (n_binf, L)

    Returns arrays:
      totals, variance, pearson_median_skew, kurtosis_excess, max_binding_offset
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
    mu2 = np.mean(centered**2, axis=1)  # sigma^2
    mu4 = np.mean(centered**4, axis=1)
    kurtosis_excess = np.zeros(counts.shape[0], dtype=np.float64)
    mask_var = mu2 > 0
    kurtosis_excess[mask_var] = mu4[mask_var] / (mu2[mask_var] ** 2) - 3.0

    # Max binding offset from sliding 5-nt window sums (center offset).
    L = counts.shape[1]
    if L < 5:
        max_binding_offset = np.zeros(counts.shape[0], dtype=np.int64)
    else:
        # sliding_window_view returns view shape (n_binf, L-4, 5)
        windows = np.lib.stride_tricks.sliding_window_view(counts, 5, axis=1)
        window_sums = windows.sum(axis=2)  # (n_binf, L-4)
        start_idx = np.argmax(window_sums, axis=1)  # 0..(L-5)
        max_binding_offset = (start_idx + 2) - window

    max_binding_offset = max_binding_offset.astype(np.int64)
    # For fully-zero loci, define max offset as 0 for interpretability.
    max_binding_offset[totals == 0] = 0

    return totals, variance, pearson_median_skew, kurtosis_excess, max_binding_offset


def smooth_metaprofile_gaussian(meta_counts: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    # Build a normalized Gaussian kernel and convolve with 'same' length output.
    radius = max(1, int(round(4.0 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()
    return np.convolve(meta_counts, kernel, mode="same")


def logistic_scale(matrix: np.ndarray) -> np.ndarray:
    # Logistic scaling after robust centering/scaling to compress large dynamic range.
    center = float(np.median(matrix))
    spread = float(np.std(matrix))
    if spread <= 0:
        spread = 1.0
    z = (matrix - center) / spread
    return 1.0 / (1.0 + np.exp(-z))


def plot_cluster_metaprofiles(
    protein_sources: list[tuple[str, Path]],
    cluster_ids_sorted: list[int],
    all_cluster_labels: np.ndarray,
    binf_windows_bed: Path,
    binf_index: dict[tuple[str, int, int], list[int]],
    window: int,
    n_binf: int,
    sigma: float,
    metaprofile_top_k: int,
    outdir: Path,
) -> None:
    offsets = np.arange(-window, window + 1, dtype=np.int64)
    cluster_masks = {cid: (all_cluster_labels == cid) for cid in cluster_ids_sorted}
    cluster_profiles: dict[int, dict[str, np.ndarray]] = {cid: {} for cid in cluster_ids_sorted}
    cluster_scores: dict[int, dict[str, float]] = {cid: {} for cid in cluster_ids_sorted}

    for protein_name, merged_xl_bed in protein_sources:
        counts = compute_counts_for_protein(
            merged_xl_bed=merged_xl_bed,
            binf_windows_bed=binf_windows_bed,
            binf_index=binf_index,
            window=window,
            n_binf=n_binf,
        )
        for cid in cluster_ids_sorted:
            mask = cluster_masks[cid]
            if not np.any(mask):
                continue
            meta_counts = counts[mask, :].mean(axis=0)
            smoothed = smooth_metaprofile_gaussian(meta_counts, sigma=sigma)
            cluster_profiles[cid][protein_name] = smoothed
            cluster_scores[cid][protein_name] = float(np.sum(smoothed))

    for cid in cluster_ids_sorted:
        profiles = cluster_profiles[cid]
        if not profiles:
            continue
        top_proteins = sorted(cluster_scores[cid], key=cluster_scores[cid].get, reverse=True)[:metaprofile_top_k]
        plt.figure(figsize=(11, 4.5))
        for protein_name in top_proteins:
            plt.plot(offsets, profiles[protein_name], label=protein_name, linewidth=2)
        plt.axvline(0, color="black", linewidth=1, alpha=0.4)
        plt.xlabel("Relative nucleotide position around inference loci (nt)")
        plt.ylabel("Mean crosslink support per region (Gaussian smoothed)")
        plt.title(f"Cluster C{cid} metaprofile")
        plt.xlim(-window, window)
        plt.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)
        plt.tight_layout(rect=[0.0, 0.0, 0.82, 1.0])
        out_path = outdir / f"metaprofile_cluster_C{cid}.png"
        plt.savefig(out_path, dpi=200)
        plt.close()
        print(f"Wrote cluster metaprofile plot to: {out_path}")


def make_tsne_from_binf_summary(
    tsv_path: Path,
    out_path: Path,
    perplexity: float,
    random_state: int,
    row_clusters: np.ndarray | None = None,
    feature_proteins: list[str] | None = None,
) -> None:
    try:
        from sklearn.manifold import TSNE
    except ImportError as e:
        raise ImportError("scikit-learn is required for --tsne. Install scikit-learn in your environment.") from e

    with open(tsv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Summary table has no header: {tsv_path}")
        if feature_proteins:
            total_cols = []
            for p in feature_proteins:
                col = f"{p}_total_overlaps"
                if col not in reader.fieldnames:
                    raise ValueError(f"Expected column {col} not found in {tsv_path}")
                total_cols.append(col)
        else:
            total_cols = [c for c in reader.fieldnames if c.endswith("_total_overlaps")]
        if not total_cols:
            raise ValueError(f"No *_total_overlaps columns found in {tsv_path}")
        rows = list(reader)

    if len(rows) < 2:
        raise ValueError("Need at least 2 rows in binf_summary.tsv for tSNE.")

    matrix = np.array([[float(r[c]) for c in total_cols] for r in rows], dtype=np.float64)
    # Keep tSNE feature scaling consistent with heatmap clustering.
    matrix = logistic_scale(matrix)
    max_perplexity = max(1.0, float(len(rows) - 1))
    used_perplexity = min(float(perplexity), max_perplexity)
    embedding = TSNE(n_components=2, perplexity=used_perplexity, random_state=random_state).fit_transform(matrix)

    plt.figure(figsize=(7, 6))
    if row_clusters is not None and len(row_clusters) == embedding.shape[0]:
        not_in_heatmap = row_clusters < 0
        low_signal = row_clusters == 0
        in_cluster = row_clusters > 0
        if np.any(not_in_heatmap):
            plt.scatter(
                embedding[not_in_heatmap, 0],
                embedding[not_in_heatmap, 1],
                s=10,
                alpha=0.25,
                c="lightgray",
                edgecolors="none",
                label="Not in heatmap filter",
            )
        if np.any(low_signal):
            plt.scatter(
                embedding[low_signal, 0],
                embedding[low_signal, 1],
                s=10,
                alpha=0.35,
                c="gray",
                edgecolors="none",
                label="C0 low-signal",
            )
        if np.any(in_cluster):
            sc = plt.scatter(
                embedding[in_cluster, 0],
                embedding[in_cluster, 1],
                s=12,
                alpha=0.85,
                c=row_clusters[in_cluster],
                cmap="tab10",
                edgecolors="none",
            )
            plt.colorbar(sc, label="Heatmap row cluster")
        if np.any(not_in_heatmap) or np.any(low_signal):
            plt.legend(loc="best", frameon=False, fontsize=8)
        if not (np.any(not_in_heatmap) or np.any(low_signal) or np.any(in_cluster)):
            plt.scatter(embedding[:, 0], embedding[:, 1], s=12, alpha=0.8, edgecolors="none")
    else:
        plt.scatter(embedding[:, 0], embedding[:, 1], s=12, alpha=0.8, edgecolors="none")
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    n_feat = len(total_cols)
    plt.title(f"tSNE ({n_feat} protein features, logistic-scaled)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Wrote tSNE plot to: {out_path} (perplexity={used_perplexity})")


def main():
    args = parse_args()
    if args.n_clusters < 1:
        raise ValueError("--n-clusters must be >= 1")
    if args.cluster_top_proteins < 1:
        raise ValueError("--cluster-top-proteins must be >= 1")
    if args.metaprofile_top_proteins < 1:
        raise ValueError("--metaprofile-top-proteins must be >= 1")
    if args.tsne and not args.table:
        raise ValueError("--tsne requires --table so binf_summary.tsv is available.")

    xldir = Path(args.xldir)
    if not xldir.exists():
        raise FileNotFoundError(f"XL directory not found: {xldir}")
    if not xldir.is_dir():
        raise NotADirectoryError(f"XL path is not a directory: {xldir}")

    binf_path = Path(args.bed)
    if not binf_path.exists():
        raise FileNotFoundError(f"Inference BED not found: {binf_path}")
    if not Path(args.genome).exists():
        raise FileNotFoundError(f"Genome sizes file not found: {args.genome}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.skip_merge:
        if args.samplesheet is not None:
            protein_sources = uniquify_names(load_samplesheet_inputs(Path(args.samplesheet), xldir))
        else:
            protein_sources = uniquify_names(discover_bed_inputs_from_xldir(xldir))
    else:
        if args.samplesheet is not None:
            print("Warning: --samplesheet is ignored unless --skip-merge is set.")
        protein_dirs: list[Path] = []
        protein_dirs.extend(list_protein_dirs(xldir))
        merged_sources: list[tuple[str, Path]] = []
        for protein_dir in protein_dirs:
            merged_dir = protein_dir.parent / "merged"
            merged_dir.mkdir(parents=True, exist_ok=True)
            merged_sources.append((protein_dir.name, merge_protein_crosslinks(protein_dir, merged_dir=merged_dir)))
        protein_sources = uniquify_names(merged_sources)

    protein_names = [name for name, _ in protein_sources]
    metaprofile_top_k = min(args.metaprofile_top_proteins, len(protein_names))

    with tempfile.TemporaryDirectory(prefix="intersect_binf_") as tmp:
        tmpdir = Path(tmp)
        binf_keys, binf_index, binf_windows_bed = load_binf_and_prepare_windows(
            binf_path=binf_path,
            window=args.window,
            genome=args.genome,
            tmpdir=tmpdir,
        )

        n_binf = len(binf_keys)
        L = 2 * args.window + 1
        offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)

        # For table output: store metrics per protein.
        totals_by_protein = {}
        variance_by_protein = {}
        skew_by_protein = {}
        kurt_by_protein = {}
        maxoff_by_protein = {}
        inspect_counts_matrix = None

        meta_profiles = {}
        meta_profile_scores = {}

        plt.figure(figsize=(11, 4.5))
        for protein_name, merged_xl_bed in protein_sources:
            counts = compute_counts_for_protein(
                merged_xl_bed=merged_xl_bed,
                binf_windows_bed=binf_windows_bed,
                binf_index=binf_index,
                window=args.window,
                n_binf=n_binf,
            )

            totals, variance, pearson_median_skew, kurtosis_excess, max_binding_offset = compute_summary_stats(
                counts=counts,
                window=args.window,
            )

            totals_by_protein[protein_name] = totals
            variance_by_protein[protein_name] = variance
            skew_by_protein[protein_name] = pearson_median_skew
            kurt_by_protein[protein_name] = kurtosis_excess
            maxoff_by_protein[protein_name] = max_binding_offset
            if args.inspect_protein == protein_name:
                inspect_counts_matrix = counts.copy()

            # Average support profile across all input binf regions.
            meta_counts = counts.mean(axis=0)
            smoothed = smooth_metaprofile_gaussian(meta_counts, sigma=args.gaussian_sigma)
            meta_profiles[protein_name] = smoothed
            # Rank proteins by total metaprofile signal and plot top-K only.
            meta_profile_scores[protein_name] = float(np.sum(smoothed))

        top_proteins = sorted(meta_profile_scores, key=meta_profile_scores.get, reverse=True)[:metaprofile_top_k]
        for protein_name in top_proteins:
            plt.plot(offsets, meta_profiles[protein_name], label=protein_name, linewidth=2)

        plt.axvline(0, color="black", linewidth=1, alpha=0.4)
        plt.xlabel("Relative nucleotide position around inference loci (nt)")
        plt.ylabel("Mean crosslink support per region (Gaussian smoothed)")
        plt.xlim(-args.window, args.window)
        plt.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)
        plt.tight_layout(rect=[0.0, 0.0, 0.82, 1.0])

        plot_path = outdir / "metaprofile.png"
        plt.savefig(plot_path, dpi=200)
        print(f"Wrote metaprofile plot to: {plot_path}")
        plt.close()

        # Heatmap: rank proteins by total signal, keep top-K for clustering/tSNE.
        protein_signal_rank = sorted(
            protein_names,
            key=lambda pn: float(np.sum(totals_by_protein[pn])),
            reverse=True,
        )
        k_top = min(args.cluster_top_proteins, len(protein_names))
        cluster_protein_names = protein_signal_rank[:k_top]
        print(
            f"Heatmap/clustering/tSNE use top {k_top} XL groups by summed total_overlaps "
            f"(of {len(protein_names)}): {', '.join(cluster_protein_names[:5])}"
            + (" ..." if k_top > 5 else "")
        )

        heatmap_matrix = np.column_stack([totals_by_protein[pn] for pn in cluster_protein_names])
        # Filter out low-support rows for stable clustering (sum across all proteins in full table).
        full_row_sums = np.column_stack([totals_by_protein[pn] for pn in protein_names]).sum(axis=1)
        heatmap_row_sums = full_row_sums
        heatmap_keep_mask = heatmap_row_sums >= 40
        heatmap_matrix_filtered = heatmap_matrix[heatmap_keep_mask, :]
        print(
            f"Global heatmap row filter (sum across all proteins >= 40): "
            f"kept {int(np.sum(heatmap_keep_mask))} / {len(heatmap_keep_mask)}"
        )
        if heatmap_matrix_filtered.shape[0] == 0:
            raise ValueError("No inference BED rows pass the global heatmap filter (sum across proteins >= 40).")

        # Logistic-scale for heatmap visualization and clustering stability.
        heatmap_matrix_scaled = logistic_scale(heatmap_matrix_filtered)
        n_filt = heatmap_matrix_scaled.shape[0]

        # Cosine distance is undefined for zero vectors; only cluster active rows.
        row_active_mask = heatmap_matrix_scaled.max(axis=1) > 0.6
        n_active = int(np.sum(row_active_mask))
        print(f"Cosine clustering active rows (max scaled support > 0.6): {n_active} / {n_filt}")
        if n_active == 0:
            raise ValueError("No active rows available for cosine clustering (all rows <= 0.6 scaled support).")
        active_matrix = heatmap_matrix_scaled[row_active_mask, :]

        n_row_clusters = min(args.n_clusters, n_active)
        if n_active > 1:
            row_distances = pdist(active_matrix, metric="cosine")
            row_linkage = linkage(row_distances, method="average")
            row_clusters = fcluster(row_linkage, t=n_row_clusters, criterion="maxclust")
        else:
            row_linkage = None
            row_clusters = np.array([1], dtype=np.int64)

        if active_matrix.shape[1] > 1:
            col_distances = pdist(active_matrix.T, metric="cosine")
            col_linkage = linkage(col_distances, method="average")
        else:
            col_linkage = None

        cluster_ids_sorted = sorted(np.unique(row_clusters))
        cluster_palette = sns.color_palette("tab10", n_colors=max(len(cluster_ids_sorted), 1))
        cluster_to_color = {cid: cluster_palette[i] for i, cid in enumerate(cluster_ids_sorted)}
        row_colors = [cluster_to_color[int(cid)] for cid in row_clusters]

        heatmap_fig = sns.clustermap(
            active_matrix,
            method="average",
            metric="cosine",
            cmap="viridis",
            row_cluster=(n_active > 1),
            col_cluster=(active_matrix.shape[1] > 1),
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            row_colors=row_colors,
            xticklabels=cluster_protein_names,
            yticklabels=False,
            figsize=(9, 10),
        )
        heatmap_fig.ax_heatmap.set_xlabel("Proteins")
        heatmap_fig.ax_heatmap.set_ylabel("Inference BED loci")
        legend_handles = [
            Patch(facecolor=cluster_to_color[cid], edgecolor="none", label=f"C{cid}") for cid in cluster_ids_sorted
        ]
        heatmap_fig.ax_heatmap.legend(
            handles=legend_handles,
            title="Row clusters",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            frameon=False,
        )
        heatmap_path = outdir / "binf_support_heatmap.png"
        heatmap_fig.savefig(heatmap_path, dpi=200)
        plt.close(heatmap_fig.fig)
        print(f"Wrote clustered heatmap to: {heatmap_path}")
        cluster_sizes = {int(cid): int(np.sum(row_clusters == cid)) for cid in np.unique(row_clusters)}
        print(f"Row cluster sizes (cosine+average): {cluster_sizes}")

        # Export per-locus cluster assignment aligned to original inference BED order.
        all_cluster_labels = np.full(n_binf, -1, dtype=np.int64)
        keep_indices = np.flatnonzero(heatmap_keep_mask)
        active_indices = keep_indices[row_active_mask]
        inactive_indices = keep_indices[~row_active_mask]
        all_cluster_labels[inactive_indices] = 0  # C0 = low-signal filtered from cosine clustering
        all_cluster_labels[active_indices] = row_clusters
        cluster_tsv_path = outdir / "binf_heatmap_clusters.tsv"
        with open(cluster_tsv_path, "w", encoding="utf-8") as fout:
            fout.write(
                "\t".join(
                    [
                        "binf_chr_start_end",
                        "chrom",
                        "start",
                        "end",
                        "row_sum_support",
                        "passes_heatmap_filter",
                        "heatmap_cluster",
                    ]
                )
                + "\n"
            )
            for i, key in enumerate(binf_keys):
                chrom, start_str, end_str = key.split("_", 2)
                passes = bool(heatmap_keep_mask[i])
                if not passes:
                    cluster_label = "NA"
                elif int(all_cluster_labels[i]) == 0:
                    cluster_label = "C0"
                else:
                    cluster_label = str(int(all_cluster_labels[i]))
                fout.write(
                    "\t".join(
                        [
                            key,
                            chrom,
                            start_str,
                            end_str,
                            f"{heatmap_row_sums[i]:.6g}",
                            "True" if passes else "False",
                            cluster_label,
                        ]
                    )
                    + "\n"
                )
        print(f"Wrote heatmap cluster assignments to: {cluster_tsv_path}")
        if args.cluster_metaprofiles:
            plot_cluster_metaprofiles(
                protein_sources=protein_sources,
                cluster_ids_sorted=cluster_ids_sorted,
                all_cluster_labels=all_cluster_labels,
                binf_windows_bed=binf_windows_bed,
                binf_index=binf_index,
                window=args.window,
                n_binf=n_binf,
                sigma=args.gaussian_sigma,
                metaprofile_top_k=metaprofile_top_k,
                outdir=outdir,
            )

        if args.inspect_protein is not None:
            if inspect_counts_matrix is None:
                available = ", ".join(protein_names)
                raise ValueError(f"Protein '{args.inspect_protein}' not found. Available proteins: {available}")
            nt_offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)
            inspect_totals = inspect_counts_matrix.sum(axis=1)
            inspect_keep_mask = inspect_totals >= 10
            inspect_counts_filtered = inspect_counts_matrix[inspect_keep_mask, :]
            inspect_counts_log = np.log1p(inspect_counts_filtered)
            print(
                f"{args.inspect_protein} inspect heatmap row filter (sum across nt >= 10): "
                f"kept {int(np.sum(inspect_keep_mask))} / {len(inspect_keep_mask)}"
            )
            protein_nt_fig = sns.clustermap(
                inspect_counts_log,
                method="average",
                metric="euclidean",
                cmap="mako",
                row_cluster=(inspect_counts_log.shape[0] > 1),
                col_cluster=False,
                xticklabels=nt_offsets,
                yticklabels=False,
                figsize=(12, 10),
            )
            protein_nt_fig.ax_heatmap.set_xlabel("Nucleotide position relative to inference locus")
            protein_nt_fig.ax_heatmap.set_ylabel("Inference BED loci")
            protein_nt_path = outdir / f"binf_{args.inspect_protein}_nt_support_heatmap.png"
            protein_nt_fig.savefig(protein_nt_path, dpi=200)
            plt.close(protein_nt_fig.fig)
            print(f"Wrote protein nucleotide heatmap to: {protein_nt_path}")

        if args.table:
            tsv_path = outdir / "binf_summary.tsv"
            header = ["binf_chr_start_end"]
            for protein_name in protein_names:
                header.extend(
                    [
                        f"{protein_name}_total_overlaps",
                        f"{protein_name}_variance",
                        f"{protein_name}_pearson_median_skew",
                        f"{protein_name}_kurtosis_excess",
                        f"{protein_name}_max_binding_offset",
                    ]
                )

            with open(tsv_path, "w", encoding="utf-8") as fout:
                fout.write("\t".join(header) + "\n")
                for i in range(n_binf):
                    row = [binf_keys[i]]
                    for protein_name in protein_names:
                        row.append(f"{totals_by_protein[protein_name][i]:.6g}")
                        row.append(f"{variance_by_protein[protein_name][i]:.6g}")
                        row.append(f"{skew_by_protein[protein_name][i]:.6g}")
                        row.append(f"{kurt_by_protein[protein_name][i]:.6g}")
                        row.append(str(int(maxoff_by_protein[protein_name][i])))
                    fout.write("\t".join(row) + "\n")

            print(f"Wrote binf summary table to: {tsv_path}")
            if args.tsne:
                tsne_path = outdir / "binf_summary_tsne.png"
                make_tsne_from_binf_summary(
                    tsv_path=tsv_path,
                    out_path=tsne_path,
                    perplexity=args.tsne_perplexity,
                    random_state=args.tsne_random_state,
                    row_clusters=all_cluster_labels,
                    feature_proteins=cluster_protein_names,
                )


if __name__ == "__main__":
    main()

