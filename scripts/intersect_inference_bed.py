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

DEFAULT_GENOME = "/camp/home/jonesm6/home/shared/genomes/hg38/hg38.genome"

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
    p.add_argument("--window", type=int, default=50, help="Half-window size in bp around each inference locus")
    p.add_argument("--gaussian-sigma", type=float, default=2.0, help="Gaussian smoothing sigma for metaprofile")
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


def main():
    args = parse_args()

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
            # Rank proteins by total metaprofile signal and plot top 10 only.
            meta_profile_scores[protein_name] = float(np.sum(smoothed))

        top_proteins = sorted(meta_profile_scores, key=meta_profile_scores.get, reverse=True)[:10]
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

        # Heatmap input matrix:
        # rows = inference loci, columns = proteins, values = total support per locus.
        heatmap_matrix = np.column_stack([totals_by_protein[pn] for pn in protein_names])
        # Filter out low-support rows for stable clustering.
        heatmap_row_sums = heatmap_matrix.sum(axis=1)
        heatmap_keep_mask = heatmap_row_sums >= 10
        heatmap_matrix_filtered = heatmap_matrix[heatmap_keep_mask, :]
        print(
            f"Global heatmap row filter (sum across proteins >= 10): "
            f"kept {int(np.sum(heatmap_keep_mask))} / {len(heatmap_keep_mask)}"
        )
        # Logistic scaling for visualization.
        heatmap_matrix_scaled = logistic_scale(heatmap_matrix_filtered)
        heatmap_fig = sns.clustermap(
            heatmap_matrix_scaled,
            method="average",
            metric="euclidean",
            cmap="viridis",
            row_cluster=(heatmap_matrix_scaled.shape[0] > 1),
            col_cluster=(len(protein_names) > 1),
            xticklabels=protein_names,
            yticklabels=False,
            figsize=(9, 10),
        )
        heatmap_fig.ax_heatmap.set_xlabel("Proteins")
        heatmap_fig.ax_heatmap.set_ylabel("Inference BED loci")
        heatmap_path = outdir / "binf_support_heatmap.png"
        heatmap_fig.savefig(heatmap_path, dpi=200)
        plt.close(heatmap_fig.fig)
        print(f"Wrote clustered heatmap to: {heatmap_path}")

        if args.inspect_protein is not None:
            if inspect_counts_matrix is None:
                available = ", ".join(protein_names)
                raise ValueError(f"Protein '{args.inspect_protein}' not found. Available proteins: {available}")
            nt_offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)
            inspect_totals = inspect_counts_matrix.sum(axis=1)
            inspect_keep_mask = inspect_totals >= 25
            inspect_counts_filtered = inspect_counts_matrix[inspect_keep_mask, :]
            print(
                f"{args.inspect_protein} inspect heatmap row filter (sum across nt >= 25): "
                f"kept {int(np.sum(inspect_keep_mask))} / {len(inspect_keep_mask)}"
            )
            protein_nt_fig = sns.clustermap(
                inspect_counts_filtered,
                method="average",
                metric="euclidean",
                cmap="mako",
                row_cluster=(inspect_counts_filtered.shape[0] > 1),
                col_cluster=(len(nt_offsets) > 1),
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


if __name__ == "__main__":
    main()

