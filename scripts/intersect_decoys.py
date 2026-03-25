import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_GENOME = "/scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38.p14-GencodeRelease44/hg38.genome"


def parse_args():
    p = argparse.ArgumentParser(description="Intersect decoys with multiple crosslink types and summarize")
    p.add_argument(
        "-x",
        "--xldir",
        required=True,
        help="Directory containing a subdirectory per protein (e.g. PRPF8/). If no subdirectories match, xldir itself is treated as one protein.",
    )
    p.add_argument(
        "-d",
        "--decoys",
        required=True,
        help="Decoy BED (must have at least 6 columns; strand in column 6).",
    )
    p.add_argument("--window", type=int, default=50, help="Half-window size in bp around each decoy")
    p.add_argument("--min-sum", type=float, default=25, help="Minimum summed signal threshold")
    p.add_argument("--genome", default=DEFAULT_GENOME, help="Genome sizes file for bedtools slop")
    p.add_argument(
        "--table",
        action="store_true",
        help="If set, write per-decoy summary TSV (always computed for all decoys; metaprofile uses filtering).",
    )
    p.add_argument(
        "-o",
        "--outdir",
        default="results",
        help="Output directory for metaprofile and (optional) TSV (default: results/ in working directory)",
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


def merge_protein_crosslinks(protein_dir: Path, merged_dir: Path) -> Path:
    """
    Merge all *genome.xl.bed files in a protein directory into a single BED6:
      chr start end . sumScore strand
    Chromosomes are normalized to chr-prefixed naming.
    """
    genome_xl = sorted(protein_dir.glob("*genome.xl.bed"))
    if not genome_xl:
        raise FileNotFoundError(f"No *genome.xl.bed found in {protein_dir}")

    merged_dir.mkdir(parents=True, exist_ok=True)
    protein_name = protein_dir.name
    out_path = merged_dir / f"{protein_name}_merged.bed"

    # Important: bedtools groupby requires sorted input by group keys.
    merge_cmd = (
        "cat *.genome.xl.bed | "
        "sort -k1,1 -k2,2n -k3,3n -k6,6 | "
        "bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | "
        "awk 'BEGIN{OFS=\"\\t\"} "
        "{chrom=$1; if (chrom !~ /^chr/) chrom=\"chr\" chrom; print chrom, $2, $3, \".\", $5, $4}' "
        f"> \"{out_path}\""
    )
    subprocess.run(merge_cmd, cwd=str(protein_dir), shell=True, check=True)
    return out_path


def load_decoys_and_prepare_windows(decoys_path: Path, window: int, genome: str, tmpdir: Path):
    """
    Writes:
      - decoys_site.bed: original decoys with extra columns containing decoy site start/end
      - decoys_windows.bed: bedtools slop output of decoys_site.bed (slopped intervals kept with extra columns)

    Returns:
      decoy_keys: list of chr_start_end strings in stable file order
      decoy_index: dict mapping (chr, start, end) -> row index
      decoys_windows_path: path to slopped decoys BED
    """
    decoys_site = tmpdir / "decoys_site.bed"
    decoys_windows = tmpdir / "decoys_windows.bed"

    decoy_keys: list[str] = []
    # Allows duplicate decoy loci (same chr/start/end can appear multiple times).
    # Maps (chr, start, end) -> list of decoy row indices that share these coordinates.
    decoy_index: dict[tuple[str, int, int], list[int]] = {}
    decoy_idx = 0

    with open(decoys_path, "r", encoding="utf-8") as fin, open(decoys_site, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                raise ValueError(f"Decoy BED has <6 columns: {line[:120]}")
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]
            # columns 4-5 are preserved (name and spliceAI score), plus we add site start/end
            key = f"{chrom}_{start}_{end}"
            decoy_keys.append(key)
            decoy_index.setdefault((chrom, start, end), []).append(decoy_idx)
            # Ensure strand stays in column 6.
            fout.write("\t".join([chrom, str(start), str(end), cols[3], cols[4], strand, str(start), str(end)]) + "\n")
            decoy_idx += 1

    subprocess.run(
        ["bedtools", "slop", "-i", str(decoys_site), "-g", genome, "-b", str(window)],
        stdout=open(decoys_windows, "w", encoding="utf-8"),
        check=True,
    )

    return decoy_keys, decoy_index, decoys_windows


def compute_counts_for_protein(
    merged_xl_bed: Path,
    decoys_windows_bed: Path,
    decoy_index: dict[tuple[str, int, int], list[int]],
    window: int,
    n_decoys: int,
):
    """
    Build a matrix counts[decoy_idx, offset_bin] where:
      offset_bin corresponds to offset in [-window..+window]
      value = summed crosslink score (merged_xl column 5) at xl_start for that offset
    Uses strand-aware overlaps (-s) and assigns offsets with decoy strand flipped so that offsets always match the decoy's 5'->3' direction.
    """
    L = 2 * window + 1
    counts = np.zeros((n_decoys, L), dtype=np.float32)

    # intersect output columns:
    # merged_xl: 6 columns
    # decoy_windows: 8 columns (chr,start,end,name,score,strand,site_start,site_end)
    # plus overlap length (-wo)
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(merged_xl_bed),
        "-b",
        str(decoys_windows_bed),
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
        decoy_chr = cols[6]
        decoy_strand = cols[11]
        decoy_site_start = int(cols[12])
        decoy_site_end = int(cols[13])
        idx_list = decoy_index.get((decoy_chr, decoy_site_start, decoy_site_end))
        if not idx_list:
            # Chromosome naming mismatch should not happen, but skip if it does.
            continue

        offset = xl_start - decoy_site_start
        if decoy_strand == "-":
            offset = -offset

        if -window <= offset <= window:
            bin_idx = offset + window
            for idx in idx_list:
                counts[idx, bin_idx] += xl_score

    _, stderr = proc.communicate()
    if proc.returncode not in (0, None):
        raise RuntimeError(f"bedtools intersect failed (code={proc.returncode}): {stderr[:500]}")

    return counts


def compute_summary_stats(counts: np.ndarray, window: int, min_sum: float):
    """
    counts shape: (n_decoys, L)

    Returns arrays:
      totals, pearson_median_skew, kurtosis_excess, max_binding_offset
    """
    totals = counts.sum(axis=1)
    mean = counts.mean(axis=1)
    median = np.median(counts, axis=1)
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
        # sliding_window_view returns view shape (n_decoys, L-4, 5)
        windows = np.lib.stride_tricks.sliding_window_view(counts, 5, axis=1)
        window_sums = windows.sum(axis=2)  # (n_decoys, L-4)
        start_idx = np.argmax(window_sums, axis=1)  # 0..(L-5)
        max_binding_offset = (start_idx + 2) - window

    max_binding_offset = max_binding_offset.astype(np.int64)
    # For fully-zero decoys, define max offset as 0 for interpretability.
    max_binding_offset[totals == 0] = 0

    return totals, pearson_median_skew, kurtosis_excess, max_binding_offset


def smooth_metaprofile(meta_counts: np.ndarray, smooth_window: int = 5) -> np.ndarray:
    kernel = np.ones(smooth_window, dtype=np.float64)
    # centered rolling sum; 'same' keeps length.
    return np.convolve(meta_counts, kernel, mode="same")


def main():
    args = parse_args()

    xldir = Path(args.xldir)
    decoys_path = Path(args.decoys)
    if not decoys_path.exists():
        raise FileNotFoundError(f"Decoys BED not found: {decoys_path}")
    if not Path(args.genome).exists():
        raise FileNotFoundError(f"Genome sizes file not found: {args.genome}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    merged_dir = xldir / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)

    protein_dirs = list_protein_dirs(xldir)
    protein_names = [p.name for p in protein_dirs]

    with tempfile.TemporaryDirectory(prefix="intersect_decoy_") as tmp:
        tmpdir = Path(tmp)
        decoy_keys, decoy_index, decoys_windows_bed = load_decoys_and_prepare_windows(
            decoys_path=decoys_path,
            window=args.window,
            genome=args.genome,
            tmpdir=tmpdir,
        )

        n_decoys = len(decoy_keys)
        L = 2 * args.window + 1
        offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)

        # For table output: store metrics per protein.
        totals_by_protein = {}
        skew_by_protein = {}
        kurt_by_protein = {}
        maxoff_by_protein = {}

        plt.figure(figsize=(10, 4.5))
        for protein_dir, protein_name in zip(protein_dirs, protein_names):
            merged_xl_bed = merge_protein_crosslinks(protein_dir, merged_dir=merged_dir)
            counts = compute_counts_for_protein(
                merged_xl_bed=merged_xl_bed,
                decoys_windows_bed=decoys_windows_bed,
                decoy_index=decoy_index,
                window=args.window,
                n_decoys=n_decoys,
            )

            totals, pearson_median_skew, kurtosis_excess, max_binding_offset = compute_summary_stats(
                counts=counts,
                window=args.window,
                min_sum=args.min_sum,
            )

            totals_by_protein[protein_name] = totals
            skew_by_protein[protein_name] = pearson_median_skew
            kurt_by_protein[protein_name] = kurtosis_excess
            maxoff_by_protein[protein_name] = max_binding_offset

            # Metaprofile uses filtering: include only decoys with totals >= min_sum.
            mask_keep = totals >= args.min_sum
            meta_counts = counts[mask_keep, :].sum(axis=0) if np.any(mask_keep) else np.zeros(L, dtype=np.float64)
            smoothed = smooth_metaprofile(meta_counts, smooth_window=5)
            plt.plot(offsets, smoothed, label=protein_name, linewidth=2)

        plt.axvline(0, color="black", linewidth=1, alpha=0.4)
        plt.xlabel("Relative nucleotide position around decoy (nt)")
        plt.ylabel("Crosslinks (smoothed sum over 5 nt)")
        plt.xlim(-args.window, args.window)
        plt.legend(frameon=False)
        plt.tight_layout()

        plot_path = outdir / "metaprofile.png"
        plt.savefig(plot_path, dpi=200)
        print(f"Wrote metaprofile plot to: {plot_path}")

        if args.table:
            tsv_path = outdir / "decoy_summary.tsv"
            header = ["decoy_chr_start_end"]
            for protein_name in protein_names:
                header.extend(
                    [
                        f"{protein_name}_total_overlaps",
                        f"{protein_name}_pearson_median_skew",
                        f"{protein_name}_kurtosis_excess",
                        f"{protein_name}_max_binding_offset",
                    ]
                )

            with open(tsv_path, "w", encoding="utf-8") as fout:
                fout.write("\t".join(header) + "\n")
                for i in range(n_decoys):
                    row = [decoy_keys[i]]
                    for protein_name in protein_names:
                        row.append(f"{totals_by_protein[protein_name][i]:.6g}")
                        row.append(f"{skew_by_protein[protein_name][i]:.6g}")
                        row.append(f"{kurt_by_protein[protein_name][i]:.6g}")
                        row.append(str(int(maxoff_by_protein[protein_name][i])))
                    fout.write("\t".join(row) + "\n")

            print(f"Wrote decoy summary table to: {tsv_path}")


if __name__ == "__main__":
    main()

