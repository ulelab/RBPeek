import argparse
import csv
import gzip
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

# Categorical hues for cluster identity, shared by the heatmap colour bar and the tSNE so a
# cluster is the same colour in both. NOT tab20, which was here before: tab20 is built as
# light/dark PAIRS, so consecutive clusters differed mainly in lightness. Measured on the
# first six slots (OKLab dE x100, adjacent pairs): tab20 scores 13.5 on normal vision against
# a floor of 15 - a fail - and tab10 scores 0.7 under simulated deuteranopia. This order
# scores 9.1 (CVD, target 8) and 20.8 (normal vision).
#
# Caveat worth keeping in mind: those are ADJACENT-pair figures. In a dense intermixed
# scatter every pair is effectively adjacent, and no six-hue set clears the all-pairs bar
# here (this one drops to 2.6 under CVD). If clusters still read as one mass, the fix is to
# facet - one panel per cluster against grey - not more hues.
CLUSTER_HUES = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # green
    "#eda100",  # amber
    "#8a63d2",  # violet
    "#00a3b5",  # teal
    "#c2477f",  # magenta
    "#6b7280",  # slate
]

DEFAULT_GENOME = "/camp/home/jonesm6/home/shared/genomes/hg38/hg38.genome"
HEATMAP_TOP_PROTEINS = 60
METAPROFILE_TOP_PROTEINS = 10
ROW_KMEANS_DEFAULT = 20

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
        "--panel-anchor",
        choices=["start", "midpoint"],
        default="start",
        help=(
            "Which point of each --xldir interval carries its score. 'start' (default) is "
            "correct for 1 nt crosslink sites, where start and midpoint are the same "
            "position. Use 'midpoint' when the panel holds PEAKS: a peak's whole score is "
            "deposited on one offset, so scoring at the start puts a peak of width w at "
            "offset -w/2 on + strand loci and +w/2 on - strand loci (the strand flip turns "
            "the genomic start into the 3' end), producing a spurious +/-w/2 doublet "
            "instead of a peak at 0."
        ),
    )
    p.add_argument(
        "--heatmap-scale",
        choices=["logistic", "percentile"],
        default="logistic",
        help=(
            "Colour scaling for the clustered locus x protein heatmap. 'logistic' (default) "
            "centres on the matrix median, which is ~0 on a sparse matrix, so every empty "
            "cell lands at exactly 0.5 and the palette below 0.5 goes unused. 'percentile' "
            "applies log1p, then scales against a high percentile of the non-zero values and "
            "clips, leaving empty cells at 0 so the full range carries signal. This also "
            "affects the column dendrogram, which is built from cosine distances on the "
            "scaled matrix."
        ),
    )
    p.add_argument(
        "--heatmap-scale-percentile",
        type=float,
        default=99.0,
        help=(
            "Percentile of non-zero values mapped to the top of the colour range under "
            "--heatmap-scale percentile (default 99). Lower it if the heatmap still looks "
            "dark, raise it if a few cells dominate."
        ),
    )
    p.add_argument(
        "--protein-select",
        choices=["total", "centrality", "enrichment"],
        default="total",
        help=(
            "How the top-K panel columns are chosen for the heatmap, clustering, tSNE and "
            "the proteins x nt heatmap. 'total' (default) ranks by summed total_overlaps, "
            "which tracks sequencing depth as much as any relationship to the inference "
            "loci. 'enrichment' (recommended) ranks by the fraction of a protein's "
            "signal-bearing loci whose strongest binding lands within --enrichment-window "
            "of the locus. 'centrality' ranks by Pearson correlation of the mean profile "
            "against a centred Gaussian; on a centred panel this ends up scoring baseline "
            "tidiness rather than position, which biases towards low-background assays - "
            "prefer 'enrichment'."
        ),
    )
    p.add_argument(
        "--enrichment-window",
        type=int,
        default=5,
        help=(
            "Half-width (nt) counted as 'on the locus' by --protein-select enrichment "
            "(default 5). A locus counts towards a protein's score when that protein's "
            "max_binding_offset there is within +/- this many nt."
        ),
    )
    p.add_argument(
        "--protein-min-loci",
        type=int,
        default=500,
        help=(
            "Minimum number of loci with any signal for a protein to be eligible under "
            "--protein-select enrichment or centrality (default 500). Without this the "
            "enrichment fraction is trivially 1.0 for a protein touching a single locus."
        ),
    )
    p.add_argument(
        "--protein-min-loci-frac",
        type=float,
        default=None,
        help=(
            "Express the eligibility gate as a FRACTION of the inference loci instead of an "
            "absolute count, and it scales with the set automatically. Overrides "
            "--protein-min-loci when given. Use this whenever two runs are meant to be "
            "compared: an absolute gate makes a smaller set harder to clear, so a protein "
            "can be ranked in one run and excluded from the other purely on set size, and "
            "exclusion is indistinguishable from absence in the output."
        ),
    )
    p.add_argument(
        "--centrality-sigma",
        type=float,
        default=5.0,
        help=(
            "Width (sigma, nt) of the Gaussian template used by --protein-select centrality "
            "(default 5). Roughly the footprint width to select for: smaller favours sharp "
            "spikes at 0, larger favours broad central humps."
        ),
    )
    p.add_argument(
        "--centrality-min-total",
        type=float,
        default=100.0,
        help=(
            "Minimum summed total_overlaps for a protein to be eligible under "
            "--protein-select centrality (default 100). Guards against a protein whose "
            "profile rests on a handful of overlaps scoring a near perfect correlation from "
            "a single spike at 0. Raise it if the selected columns look sparse."
        ),
    )
    p.add_argument(
        "--exclude-groups",
        default=None,
        help=(
            "Comma-separated substrings; any panel column whose group name contains one is "
            "dropped before anything is computed. E.g. --exclude-groups PARCLIP. Note this "
            "removes those columns from the row-support sums behind --heatmap-min-support "
            "as well, so that threshold may need retuning alongside it."
        ),
    )
    p.add_argument(
        "--heatmap-min-support",
        type=float,
        default=50.0,
        help=(
            "Keep inference loci whose summed support across ALL xl groups is at least this "
            "(default 50, the previously hardcoded value). This scales with the number of "
            "panel columns, not with anything intrinsic to a locus, so a wide panel makes it "
            "permissive: 87%% of rows passed on a 297-column run, which saturates the heatmap. "
            "Raise it when the heatmap looks uniformly dark. Prefer "
            "--heatmap-min-support-percentile whenever two runs are meant to be compared."
        ),
    )
    p.add_argument(
        "--heatmap-min-support-percentile",
        type=float,
        default=None,
        help=(
            "Size-matched alternative to --heatmap-min-support, which it overrides when given. "
            "Drops the bottom P%% of loci by summed support instead of cutting at an absolute "
            "value, so two runs over different locus sets are filtered equally hard. Row "
            "support has no intrinsic scale - it is the summed SCORE of every overlapping "
            "panel interval across every column - so an absolute threshold means different "
            "things in different runs: at 200, the THRAP3 exonic run lost 11.3%% of loci and "
            "the intronic run 43.5%%, because median support differs 5x between them. Loci "
            "with zero support are ALWAYS dropped, whatever the percentile resolves to; "
            "without that, a set where over P%% of loci are empty (18.6%% of the THRAP3 "
            "intronic loci) yields a threshold of 0 and plots blank heatmap rows."
        ),
    )
    p.add_argument(
        "--nt-heatmap-window",
        type=int,
        default=None,
        help=(
            "Crop the -i/--inspect-protein nucleotide heatmap to +/- this many nt for DISPLAY "
            "only (default: "
            "show the full --window). Purely a plotting crop: --window still governs what is "
            "counted, so every statistic, the metaprofile, the enrichment ranking and the "
            "summary table are unchanged. Use it to zoom in on the centre without narrowing "
            "the analysis, e.g. --window 100 --nt-heatmap-window 50."
        ),
    )
    p.add_argument(
        "--no-clustering",
        action="store_true",
        help=(
            "Skip k-means row clustering entirely. The loci x protein heatmap is still "
            "written, with rows ordered by total support instead of grouped by cluster, and "
            "so is the tSNE, as a single-colour scatter. Suppresses binf_heatmap_clusters.tsv "
            "and is incompatible with --cluster-metaprofiles. Use it when the question is "
            "'which panel columns are enriched here', not 'what kinds of locus are there' - "
            "the clusters are only worth computing if something downstream consumes them."
        ),
    )
    p.add_argument(
        "--cluster-metaprofiles",
        action="store_true",
        help="If set, write one metaprofile plot per heatmap row cluster.",
    )
    p.add_argument(
        "--n-clusters",
        type=int,
        default=ROW_KMEANS_DEFAULT,
        help=(
            "Number of k-means clusters on the binarized (presence/absence) heatmap row matrix "
            f"(default: {ROW_KMEANS_DEFAULT}). Capped at the number of filtered rows."
        ),
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
            "Plot the top K proteins in the global and per-cluster metaprofiles "
            f"(default: {METAPROFILE_TOP_PROTEINS}). These are the first K of the SAME "
            "ranking that picks the heatmap's columns (--protein-select), so the metaprofile "
            "is always a subset of the heatmap and the two figures agree on which proteins "
            "matter. Capped at --cluster-top-proteins."
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


def _chrom_style(path: Path, sample: int = 2000):
    """
    True if this BED is chr-prefixed, False if Ensembl-style, None if it has no data rows.

    Decided by MAJORITY over the first `sample` data rows, not by the first row alone.
    Sampling one line is wrong for a very common case: `sort -k1,1` is ASCII, so uppercase
    scaffold names (GL000009.2, KI270302.1) sort BEFORE "chr", and a perfectly chr-prefixed
    panel file whose first row is a scaffold looks Ensembl-style. That misread 278 of 302
    columns as needing a rewrite on one run - harmless, since the rewrite only ever ADDS a
    prefix to names lacking one, but 278 pointless file copies per run.
    """
    n_chr = n_other = 0
    with _open_text_auto(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
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

    This exists because the failure it prevents is SILENT. bedtools finds no overlap between
    "1" and "chr1", so a mismatched panel column yields all zeros, the run completes happily,
    and the column simply looks like an RBP that binds nothing. Five columns were lost that
    way before this check existed.

    Only mismatched files are rewritten, into tmpdir, so a correctly-named panel costs one
    line read per file and nothing else.
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
                if not line.strip() or line.startswith(("#", "track", "browser")):
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
    panel_anchor: str = "start",
):
    """
    Build a matrix counts[binf_idx, offset_bin] where:
      offset_bin corresponds to offset in [-window..+window]
      value = summed crosslink score (merged_xl column 5) at the interval's anchor for that offset
    Uses strand-aware overlaps (-s) and assigns offsets with strand flipped so that offsets always match the inference BED's 5'->3' direction.

    panel_anchor selects which point of each panel interval carries its score. The whole
    score lands on ONE offset either way, so the choice matters as soon as intervals are
    wider than 1 nt: with "start", a width-w peak centred on a locus scores at -w/2 on
    + strand loci, and the strand flip below sends the same peak to +w/2 on - strand loci
    (genomic start is the 3' end there), splitting a real central signal into a spurious
    +/-w/2 doublet. "midpoint" collapses both back onto 0. For true 1 nt crosslink sites
    the two options are identical, since start == midpoint when end == start + 1.
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
        if panel_anchor == "midpoint":
            xl_start = (xl_start + int(cols[2])) // 2
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


def centrality_rank(
    protein_names: list[str],
    meta_profiles: dict[str, np.ndarray],
    totals_by_protein: dict[str, np.ndarray],
    window: int,
    sigma: float,
    min_total: float,
) -> tuple[list[str], dict[str, float], list[str]]:
    """
    Rank proteins by how well their mean support profile matches a Gaussian centred on the
    inference locus, i.e. by how specifically concentrated their binding is at offset 0.

    Ranking by summed total_overlaps instead (the default) tracks how deeply a dataset was
    sequenced and how many peaks it has genome-wide, which is why peak-rich panels dominate
    the top-K regardless of whether their binding has anything to do with the inference loci.
    Pearson r against a fixed template is scale-free, so a shallow dataset with a sharp
    central profile can outrank a deep one with a flat profile.

    Proteins under min_total are excluded first. Their profiles are built from very few
    overlaps, so a single locus contributing one spike at 0 would otherwise score a near
    perfect correlation. Flat profiles are excluded too, since Pearson r is undefined when
    either vector has zero variance.

    Returns (ranked_names, scores, excluded_names).
    """
    offsets = np.arange(-window, window + 1, dtype=np.float64)
    template = np.exp(-0.5 * (offsets / float(sigma)) ** 2)
    template = template - template.mean()
    template_norm = float(np.sqrt(np.sum(template**2)))

    scores: dict[str, float] = {}
    excluded: list[str] = []
    for pn in protein_names:
        if float(np.sum(totals_by_protein[pn])) < min_total:
            excluded.append(pn)
            continue
        prof = np.asarray(meta_profiles[pn], dtype=np.float64)
        centred = prof - prof.mean()
        denom = float(np.sqrt(np.sum(centred**2))) * template_norm
        if denom <= 0:
            excluded.append(pn)
            continue
        scores[pn] = float(np.dot(centred, template) / denom)

    ranked = sorted(scores, key=scores.get, reverse=True)
    return ranked, scores, excluded


def enrichment_rank(
    protein_names: list[str],
    totals_by_protein: dict[str, np.ndarray],
    maxoff_by_protein: dict[str, np.ndarray],
    win: int,
    min_total: float,
    min_loci: int,
) -> tuple[list[str], dict[str, float], list[str]]:
    """
    Rank proteins by the fraction of their signal-bearing loci whose strongest binding falls
    within +/-win nt of the inference locus.

    This is a PER-LOCUS measure, which is what makes it behave differently from
    centrality_rank's correlation against the mean profile. Once the panel is centred with
    --panel-anchor midpoint, nearly every mean profile peaks at 0, so correlation against a
    centred template stops discriminating on position and starts discriminating on how tidy
    a dataset's baseline is. That favours assays with low background: on the THRAP3 panel it
    put PAR-CLIP at 2.2x its share of the panel and dropped BCLAF1, THRAP3's known complex
    partner, out of the top 20 entirely. Asking instead "at what fraction of loci does this
    protein actually peak on the locus" is insensitive to baseline texture. On the same data
    it returned 19 eCLIP + 1 iCLIP, matching panel composition, with both HNRNPC datasets
    first and second and BCLAF1 back inside the top 20.

    Both gates matter. Without min_loci the metric is trivially 1.0 for a protein whose
    signal touches a single locus: the raw top of this ranking was a PAR-CLIP column with
    one locus and a summed score of 5.

    Returns (ranked_names, scores, excluded_names).
    """
    scores: dict[str, float] = {}
    excluded: list[str] = []
    for pn in protein_names:
        totals = np.asarray(totals_by_protein[pn])
        has_signal = totals > 0
        n_loci = int(has_signal.sum())
        if float(totals.sum()) < min_total or n_loci < min_loci:
            excluded.append(pn)
            continue
        offs = np.asarray(maxoff_by_protein[pn])[has_signal]
        scores[pn] = float(np.mean(np.abs(offs) <= win))

    ranked = sorted(scores, key=scores.get, reverse=True)
    return ranked, scores, excluded


def logistic_scale(matrix: np.ndarray) -> np.ndarray:
    # Logistic scaling after robust centering/scaling to compress large dynamic range.
    #
    # Note this centres on the matrix median, which is ~0 whenever the matrix is sparse.
    # Every empty cell then maps to exactly 0.5, so half the palette is never used and the
    # colourbar starts at 0.5 rather than 0. See percentile_scale for the alternative.
    center = float(np.median(matrix))
    spread = float(np.std(matrix))
    if spread <= 0:
        spread = 1.0
    z = (matrix - center) / spread
    return 1.0 / (1.0 + np.exp(-z))


def percentile_scale(matrix: np.ndarray, pct: float) -> np.ndarray:
    """
    log1p the matrix, then scale to [0, 1] against a high percentile of the NON-ZERO values.

    Empty cells stay at 0 rather than being pushed to mid-palette, so the whole colour range
    carries signal. The percentile is taken over non-zero values because on a sparse matrix
    a percentile over all cells is itself 0, which collapses the scale.

    The log1p step is what makes this usable rather than merely correct. Support counts are
    heavy-tailed: on a representative matrix the non-zero median was 10 against a maximum of
    333, so dividing raw values by the 99th percentile rendered the median cell at 0.11 and
    left the heatmap darker than the logistic scaling it replaced. After log1p the same
    settings put the median cell at 0.53, with ~59% of non-zero cells above mid-palette.

    Clipping at pct rather than the maximum keeps a handful of very large cells from
    compressing everything else; cells above the percentile saturate at 1.0.
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


METAPROFILE_YLABEL = "Mean peak support across all loci (Gaussian smoothed)"


def render_metaprofile(
    offsets: np.ndarray,
    profiles: dict,
    order: list,
    grand_totals: dict,
    n_loci: int,
    window: int,
    out_path: Path,
    title: str,
) -> None:
    """
    Draw one metaprofile panel.

    The left axis is mean support per locus, i.e. counts.mean(axis=0). Every locus is in the
    denominator, including the ones where the protein has no signal at all, so a curve is
    diluted by non-binding loci rather than describing the sites it does bind.

    The right axis is that multiplied by n_loci, which recovers the summed panel score.
    Because n_loci is a single constant the two axes are the same curve at two scales and
    agree pixel for pixel - the only normalisation that can honestly share an axis. Anything
    per-protein (dividing each curve by its own maximum, say) reorders the curves and needs
    its own panel.
    """
    # Height grows with the legend so 20 curves stay as legible as 10.
    fig_h = max(4.8, 0.30 * len(order) + 2.2)
    fig, ax = plt.subplots(figsize=(13.5, fig_h))
    # The default prop cycle is 10 colours, so a 20-curve plot silently draws two proteins in
    # the same blue. Vary LINESTYLE every time the colours wrap instead of reaching for a
    # 20-hue palette: tab20 does not clear the normal-vision separation floor, and a second
    # non-colour channel stays readable under any colour-vision deficiency.
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    linestyles = ["-", "--", ":", "-."]
    for i, protein_name in enumerate(order):
        total = float(grand_totals.get(protein_name, 0.0))
        ax.plot(
            offsets,
            profiles[protein_name],
            label=f"{protein_name}  (sum={total:,.0f})",
            color=palette[i % len(palette)],
            linestyle=linestyles[(i // len(palette)) % len(linestyles)],
            linewidth=2,
        )
    ax.axvline(0, color="black", linewidth=1, alpha=0.4)
    ax.set_xlabel("Relative nucleotide position around inference loci (nt)")
    ax.set_ylabel(METAPROFILE_YLABEL)
    ax.set_xlim(-window, window)
    # Left-aligned: sample names make these titles long, and a centred one overruns the
    # figure once tight_layout has shrunk the axes to make room for the legend.
    ax.set_title(title, loc="left", fontsize=10)

    ax2 = ax.twinx()
    lo, hi = ax.get_ylim()
    ax2.set_ylim(lo * n_loci, hi * n_loci)
    ax2.set_ylabel(f"Total peak support (sum over {n_loci:,} loci)")

    # Fixed FIGURE-fraction layout rather than tight_layout + an axes-relative legend. Three
    # things compete for the strip to the right of the plot: ax2's tick labels, ax2's y-label,
    # and the legend. Anchoring the legend in axes coordinates puts it a fixed multiple of the
    # AXES width away, which is blind to how wide the tick labels turned out - six-digit
    # totals pushed the y-label straight under the legend text. Figure coordinates make the
    # gap absolute, so it holds whatever the numbers are.
    fig.subplots_adjust(left=0.06, right=0.50, top=1.0 - 0.5 / fig_h, bottom=0.72 / fig_h)
    fig.legend(
        *ax.get_legend_handles_labels(),
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.60, 0.5),
        bbox_transform=fig.transFigure,
        borderaxespad=0.0,
    )
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


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
    panel_anchor: str = "start",
) -> None:
    offsets = np.arange(-window, window + 1, dtype=np.int64)
    cluster_masks = {cid: (all_cluster_labels == cid) for cid in cluster_ids_sorted}
    cluster_profiles: dict[int, dict[str, np.ndarray]] = {cid: {} for cid in cluster_ids_sorted}
    cluster_scores: dict[int, dict[str, float]] = {cid: {} for cid in cluster_ids_sorted}
    cluster_totals: dict[int, dict[str, float]] = {cid: {} for cid in cluster_ids_sorted}

    for protein_name, merged_xl_bed in protein_sources:
        counts = compute_counts_for_protein(
            merged_xl_bed=merged_xl_bed,
            binf_windows_bed=binf_windows_bed,
            binf_index=binf_index,
            window=window,
            n_binf=n_binf,
            panel_anchor=panel_anchor,
        )
        for cid in cluster_ids_sorted:
            mask = cluster_masks[cid]
            if not np.any(mask):
                continue
            meta_counts = counts[mask, :].mean(axis=0)
            smoothed = smooth_metaprofile_gaussian(meta_counts, sigma=sigma)
            cluster_profiles[cid][protein_name] = smoothed
            cluster_scores[cid][protein_name] = float(np.sum(smoothed))
            cluster_totals[cid][protein_name] = float(counts[mask, :].sum())

    for cid in cluster_ids_sorted:
        profiles = cluster_profiles[cid]
        if not profiles:
            continue
        top_proteins = sorted(cluster_scores[cid], key=cluster_scores[cid].get, reverse=True)[:metaprofile_top_k]
        n_loci = int(np.sum(cluster_masks[cid]))
        out_path = outdir / f"metaprofile_cluster_C{cid}.png"
        render_metaprofile(
            offsets=offsets,
            profiles=profiles,
            order=top_proteins,
            grand_totals=cluster_totals[cid],
            n_loci=n_loci,
            window=window,
            out_path=out_path,
            title=f"Cluster C{cid} metaprofile (n = {n_loci:,} loci)",
        )
        print(f"Wrote cluster metaprofile plot to: {out_path}")


def make_tsne_from_binf_summary(
    tsv_path: Path,
    out_path: Path,
    perplexity: float,
    random_state: int,
    row_clusters: np.ndarray | None = None,
    feature_proteins: list[str] | None = None,
    cluster_to_color: dict | None = None,
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
        if np.any(in_cluster):
            # Cluster id is CATEGORICAL. It used to be handed to plt.scatter as a numeric
            # array with cmap="tab20" and a colorbar, which renders a 20-hue qualitative map
            # as a continuous ramp: the reader gets a rainbow scale bar for what are six
            # discrete groups, and the hues do not match the heatmap's own cluster colours.
            # Draw one series per cluster instead, in the SAME colours the heatmap uses, so
            # C3 in this plot is C3 over there.
            for cid in sorted(np.unique(row_clusters[in_cluster])):
                m = row_clusters == cid
                colour = (cluster_to_color or {}).get(int(cid))
                plt.scatter(
                    embedding[m, 0],
                    embedding[m, 1],
                    s=12,
                    alpha=0.85,
                    color=colour,
                    edgecolors="none",
                    label=f"C{int(cid)}",
                )
        if np.any(not_in_heatmap) or np.any(in_cluster):
            plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
                       fontsize=8, title="Locus cluster")
        if not (np.any(not_in_heatmap) or np.any(in_cluster)):
            plt.scatter(embedding[:, 0], embedding[:, 1], s=12, alpha=0.8, edgecolors="none")
    else:
        plt.scatter(embedding[:, 0], embedding[:, 1], s=12, alpha=0.8, edgecolors="none")
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    n_feat = len(total_cols)
    plt.title(f"tSNE ({n_feat} protein features, logistic-scaled)")
    plt.tight_layout(rect=[0.0, 0.0, 0.85, 1.0])
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
    if args.no_clustering and args.cluster_metaprofiles:
        raise ValueError(
            "--cluster-metaprofiles needs the row clusters that --no-clustering suppresses; "
            "pass one or the other."
        )
    if args.heatmap_min_support_percentile is not None and not 0 < args.heatmap_min_support_percentile < 100:
        raise ValueError("--heatmap-min-support-percentile must be between 0 and 100 (exclusive)")

    xldir = Path(args.xldir)
    if not xldir.exists():
        raise FileNotFoundError(f"XL directory not found: {xldir}")
    if not xldir.is_dir():
        raise NotADirectoryError(f"XL path is not a directory: {xldir}")

    binf_path = Path(args.bed)
    if not binf_path.exists():
        raise FileNotFoundError(f"Inference BED not found: {binf_path}")
    if args.nt_heatmap_window is not None:
        if args.nt_heatmap_window < 1:
            raise ValueError("--nt-heatmap-window must be >= 1")
        if args.nt_heatmap_window > args.window:
            raise ValueError(
                f"--nt-heatmap-window {args.nt_heatmap_window} exceeds --window {args.window}; "
                "it is a display crop, so it cannot show more than was counted."
            )
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

    if args.exclude_groups:
        patterns = [p.strip() for p in args.exclude_groups.split(",") if p.strip()]
        before = len(protein_sources)
        dropped_names = [n for n, _ in protein_sources if any(p in n for p in patterns)]
        protein_sources = [(n, p) for n, p in protein_sources if not any(pat in n for pat in patterns)]
        if not protein_sources:
            raise ValueError(f"--exclude-groups {args.exclude_groups!r} removed every panel column.")
        print(
            f"Excluded {before - len(protein_sources)} of {before} panel columns matching "
            f"{patterns}: {', '.join(dropped_names[:5])}" + (" ..." if len(dropped_names) > 5 else "")
        )

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

        protein_sources = harmonise_panel_chroms(protein_sources, binf_path, tmpdir)
        protein_names = [name for name, _ in protein_sources]

        n_binf = len(binf_keys)
        offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)

        # For table output: store metrics per protein.
        totals_by_protein = {}
        variance_by_protein = {}
        skew_by_protein = {}
        kurt_by_protein = {}
        maxoff_by_protein = {}
        inspect_counts_matrix = None

        meta_profiles = {}

        for protein_name, merged_xl_bed in protein_sources:
            counts = compute_counts_for_protein(
                merged_xl_bed=merged_xl_bed,
                binf_windows_bed=binf_windows_bed,
                binf_index=binf_index,
                window=args.window,
                n_binf=n_binf,
                panel_anchor=args.panel_anchor,
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

        # The metaprofile is plotted below, AFTER the protein ranking, because it now draws
        # its top-K from that same ranking rather than from depth.

        # Resolve the eligibility gate. A fractional gate is the safe choice when two runs
        # will be compared, because an absolute one is harsher on the smaller locus set.
        min_loci_eff = args.protein_min_loci
        if args.protein_min_loci_frac is not None:
            if not 0 < args.protein_min_loci_frac <= 1:
                raise ValueError("--protein-min-loci-frac must be in (0, 1]")
            min_loci_eff = int(math.ceil(args.protein_min_loci_frac * n_binf))
            print(
                f"Eligibility gate: --protein-min-loci-frac {args.protein_min_loci_frac:g} "
                f"x {n_binf} loci = {min_loci_eff} (overrides --protein-min-loci "
                f"{args.protein_min_loci})"
            )

        # Heatmap: rank proteins, keep top-K for clustering/tSNE.
        if args.protein_select == "enrichment":
            protein_signal_rank, centrality_scores, centrality_excluded = enrichment_rank(
                protein_names=protein_names,
                totals_by_protein=totals_by_protein,
                maxoff_by_protein=maxoff_by_protein,
                win=args.enrichment_window,
                min_total=args.centrality_min_total,
                min_loci=min_loci_eff,
            )
            if not protein_signal_rank:
                raise ValueError(
                    f"No protein passed --centrality-min-total {args.centrality_min_total:g} "
                    f"and a minimum of {min_loci_eff} loci with signal. Lower them, or use "
                    "--protein-select total."
                )
            if centrality_excluded:
                print(
                    f"Enrichment ranking: {len(centrality_excluded)} of {len(protein_names)} "
                    f"proteins excluded (summed total_overlaps < {args.centrality_min_total:g} "
                    f"or fewer than {min_loci_eff} loci with signal)"
                )
            rank_desc = f"fraction of loci peaking within +/-{args.enrichment_window}nt"
        elif args.protein_select == "centrality":
            protein_signal_rank, centrality_scores, centrality_excluded = centrality_rank(
                protein_names=protein_names,
                meta_profiles=meta_profiles,
                totals_by_protein=totals_by_protein,
                window=args.window,
                sigma=args.centrality_sigma,
                min_total=args.centrality_min_total,
            )
            if not protein_signal_rank:
                raise ValueError(
                    f"No protein passed --centrality-min-total {args.centrality_min_total:g}. "
                    "Lower it, or use --protein-select total."
                )
            if centrality_excluded:
                print(
                    f"Centrality ranking: {len(centrality_excluded)} of {len(protein_names)} "
                    f"proteins excluded (summed total_overlaps < {args.centrality_min_total:g}, "
                    "or flat profile)"
                )
            rank_desc = f"Pearson r vs a sigma={args.centrality_sigma:g} Gaussian centred on the locus"
        else:
            protein_signal_rank = sorted(
                protein_names,
                key=lambda pn: float(np.sum(totals_by_protein[pn])),
                reverse=True,
            )
            centrality_scores = None
            rank_desc = "summed total_overlaps"

        k_top = min(args.cluster_top_proteins, len(protein_signal_rank))
        cluster_protein_names = protein_signal_rank[:k_top]
        print(
            f"Heatmap/clustering/tSNE use top {k_top} XL groups by {rank_desc} "
            f"(of {len(protein_names)}): {', '.join(cluster_protein_names[:5])}"
            + (" ..." if k_top > 5 else "")
        )
        score_label = {
            "enrichment": "frac centred",
            "centrality": "centrality r",
            "total": "mean peak support",
        }[args.protein_select]
        if centrality_scores is not None:
            # Print the selected columns with both scores, so a high-r/low-signal column is
            # obvious rather than silently shaping the heatmap.
            print(f"  selected columns ({score_label}, loci with signal, summed total_overlaps):")
            for pn in cluster_protein_names:
                totals = np.asarray(totals_by_protein[pn])
                print(
                    f"    {pn:40} {score_label}={centrality_scores[pn]:+.3f}  "
                    f"loci={int((totals > 0).sum()):>6}  total={float(totals.sum()):,.0f}"
                )

        # The metaprofile draws from the SAME ranking that picks the heatmap columns, so the
        # two figures never disagree about which proteins matter. Previously it ranked by
        # sum(smoothed profile), i.e. purely by depth, while the heatmap ranked by enrichment -
        # which is why a separate proteins x nt heatmap had to exist to show the heatmap's
        # proteins as profiles.
        meta_top_k = min(metaprofile_top_k, len(cluster_protein_names))
        grand_totals = {pn: float(np.sum(totals_by_protein[pn])) for pn in protein_names}
        plot_path = outdir / "metaprofile.png"
        render_metaprofile(
            offsets=offsets,
            profiles=meta_profiles,
            order=cluster_protein_names[:meta_top_k],
            grand_totals=grand_totals,
            n_loci=n_binf,
            window=args.window,
            out_path=plot_path,
            title=f"Top {meta_top_k} of {len(protein_names)} panel samples by {rank_desc}"
                  f"   |   n = {n_binf:,} loci",
        )
        print(f"Wrote metaprofile plot to: {plot_path}")

        # Ranked panel samples. This is the table that carries BOTH rankings side by side -
        # mean support (depth around the locus) and frac_centred (how central that support is) -
        # which are the two competing answers to "which RBP co-occupies these sites".
        ranking_path = outdir / "protein_ranking.tsv"
        by_mean = sorted(protein_names, key=lambda pn: grand_totals[pn], reverse=True)
        enr_rank = {pn: i + 1 for i, pn in enumerate(protein_signal_rank)}
        selected = set(cluster_protein_names)
        with open(ranking_path, "w", encoding="utf-8") as fout:
            fout.write(
                "\t".join(
                    [
                        "rank",
                        "sample",
                        "mean_peak_support",
                        "total_peak_support",
                        "loci_with_signal",
                        "frac_loci_with_signal",
                        f"{args.protein_select}_score",
                        f"{args.protein_select}_rank",
                        "eligible",
                        "selected_for_heatmap",
                    ]
                )
                + "\n"
            )
            for i, pn in enumerate(by_mean, 1):
                totals = np.asarray(totals_by_protein[pn])
                n_sig = int((totals > 0).sum())
                # Under --protein-select total the ranking score IS the summed support, so
                # report that rather than an NA column that duplicates total_peak_support.
                if centrality_scores is None:
                    score = f"{grand_totals[pn]:.6g}"
                elif pn in centrality_scores:
                    score = f"{centrality_scores[pn]:.6g}"
                else:
                    score = "NA"
                fout.write(
                    "\t".join(
                        [
                            str(i),
                            pn,
                            f"{grand_totals[pn] / n_binf:.6g}",
                            f"{grand_totals[pn]:.6g}",
                            str(n_sig),
                            f"{n_sig / n_binf:.6g}",
                            score,
                            str(enr_rank.get(pn, "NA")),
                            "True" if pn in enr_rank else "False",
                            "True" if pn in selected else "False",
                        ]
                    )
                    + "\n"
                )
        print(f"Wrote ranked panel samples to: {ranking_path} (sorted by mean peak support)")

        heatmap_matrix = np.column_stack([totals_by_protein[pn] for pn in cluster_protein_names])
        # Filter out low-support rows for stable clustering (sum across all proteins in full table).
        full_row_sums = np.column_stack([totals_by_protein[pn] for pn in protein_names]).sum(axis=1)
        heatmap_row_sums = full_row_sums
        n_zero = int(np.sum(heatmap_row_sums <= 0))
        if args.heatmap_min_support_percentile is not None:
            # Percentile of every row, then drop empty rows regardless. The zero clause is
            # load-bearing: where more than P% of loci have no support at all the percentile
            # resolves to 0 and a bare threshold would keep those loci as blank heatmap rows.
            support_threshold = float(np.percentile(heatmap_row_sums, args.heatmap_min_support_percentile))
            heatmap_keep_mask = (heatmap_row_sums > 0) & (heatmap_row_sums >= support_threshold)
            filter_desc = (
                f"bottom {args.heatmap_min_support_percentile:g}% by summed support "
                f"(threshold {support_threshold:.6g}) plus all empty rows"
            )
        else:
            support_threshold = float(args.heatmap_min_support)
            heatmap_keep_mask = heatmap_row_sums >= support_threshold
            filter_desc = f"sum across all proteins >= {support_threshold:g}"
        heatmap_matrix_filtered = heatmap_matrix[heatmap_keep_mask, :]
        n_kept = int(np.sum(heatmap_keep_mask))
        n_rows = len(heatmap_keep_mask)
        print(
            f"Global heatmap row filter ({filter_desc}): kept {n_kept} / {n_rows} "
            f"({100 * n_kept / n_rows:.1f}%), dropped {n_rows - n_kept} "
            f"({100 * (n_rows - n_kept) / n_rows:.1f}%), of which {n_zero} had zero support "
            f"({100 * n_zero / n_rows:.1f}% of all loci)"
        )
        if heatmap_matrix_filtered.shape[0] == 0:
            raise ValueError(
                f"No inference BED rows pass the global heatmap filter ({filter_desc}). "
                "Lower --heatmap-min-support or --heatmap-min-support-percentile."
            )

        # Logistic-scale for heatmap display; k-means clustering uses binarized presence/absence.
        if args.heatmap_scale == "percentile":
            heatmap_matrix_scaled = percentile_scale(heatmap_matrix_filtered, args.heatmap_scale_percentile)
            frac_saturated = float(np.mean(heatmap_matrix_scaled >= 1.0))
            print(
                f"Heatmap colour scale: percentile (clip at {args.heatmap_scale_percentile:g}th of "
                f"non-zero), {100 * float(np.mean(heatmap_matrix_filtered == 0)):.1f}% of cells empty, "
                f"{100 * frac_saturated:.2f}% saturated at 1.0"
            )
        else:
            heatmap_matrix_scaled = logistic_scale(heatmap_matrix_filtered)
        n_filt = heatmap_matrix_scaled.shape[0]
        row_totals_filtered = heatmap_row_sums[heatmap_keep_mask]

        if args.no_clustering:
            row_clusters = None
            row_clusters_sorted = None
            cluster_ids_sorted = []
            cluster_to_color = None
            row_colors = None
            k_fit = 0
            # No cluster key to sort within, so support alone orders the locus axis.
            sort_idx = np.argsort(-row_totals_filtered.astype(np.float64), kind="stable")
            heatmap_display = heatmap_matrix_scaled[sort_idx]
        else:
            binary_rows = (heatmap_matrix_filtered > 0.0).astype(np.float64)

            try:
                from sklearn.cluster import KMeans
            except ImportError as e:
                raise ImportError("scikit-learn is required for k-means row clustering. Install scikit-learn.") from e

            k_fit = min(args.n_clusters, n_filt)
            if k_fit < 1:
                raise ValueError("No rows available for k-means clustering.")
            km = KMeans(n_clusters=k_fit, random_state=args.tsne_random_state, n_init="auto")
            km_labels = km.fit_predict(binary_rows).astype(np.int64)
            # sklearn uses 0..k-1; use 1..k for downstream labels
            row_clusters = km_labels + 1

            # Row order: by cluster id (asc), then by total support (desc) within cluster.
            sort_idx = np.lexsort((-row_totals_filtered.astype(np.float64), row_clusters))
            heatmap_display = heatmap_matrix_scaled[sort_idx]
            row_clusters_sorted = row_clusters[sort_idx]

        # Column linkage only: cosine distance between proteins across rows (order-independent).
        if heatmap_matrix_scaled.shape[1] > 1:
            col_distances = pdist(heatmap_matrix_scaled.T, metric="cosine")
            col_linkage = linkage(col_distances, method="average")
        else:
            col_linkage = None

        if not args.no_clustering:
            cluster_ids_sorted = sorted(np.unique(row_clusters))
            n_colors = max(len(cluster_ids_sorted), 1)
            if n_colors <= len(CLUSTER_HUES):
                cluster_palette = CLUSTER_HUES[:n_colors]
            else:
                cluster_palette = sns.color_palette("husl", n_colors=n_colors)
            cluster_to_color = {cid: cluster_palette[i] for i, cid in enumerate(cluster_ids_sorted)}
            row_colors = [cluster_to_color[int(cid)] for cid in row_clusters_sorted]

        if args.heatmap_scale == "percentile":
            cbar_label = f"log1p support\n({args.heatmap_scale_percentile:g}th pct clip)"
            # Pin the range so empty cells read as the bottom of the palette rather than
            # being stretched by whatever the observed minimum happens to be.
            scale_kwargs = {"vmin": 0.0, "vmax": 1.0}
        else:
            cbar_label = "support\n(logistic; empty = 0.5)"
            scale_kwargs = {}
        # Transposed: proteins on the ROWS so their names read horizontally on the left,
        # instead of ~20 vertical labels along the bottom. Loci go on the x-axis, where they
        # are unlabelled anyway (there are thousands), and their cluster colour bar moves
        # from row_colors to col_colors. col_linkage was cosine distance BETWEEN PROTEINS, so
        # after the transpose it is the row linkage; loci stay in their pre-sorted
        # cluster order, so the locus axis is not re-clustered.
        n_prot = heatmap_display.shape[1]
        heatmap_fig = sns.clustermap(
            heatmap_display.T,
            row_cluster=(n_prot > 1),
            row_linkage=col_linkage,
            col_cluster=False,
            col_colors=row_colors,
            cmap="viridis",
            # Rows are ordered by the protein dendrogram, which groups by co-occurrence and
            # therefore puts the deepest columns together - not by selection rank. Printing
            # the rank makes the two orderings comparable at a glance, and makes it obvious
            # when a visually dominant row is one the ranking placed near the cut.
            yticklabels=[f"{pn}  [{i}]" for i, pn in enumerate(cluster_protein_names, 1)],
            xticklabels=False,
            figsize=(11, max(5.0, 0.34 * n_prot + 2.0)),
            cbar_kws={"label": cbar_label},
            **scale_kwargs,
        )
        heatmap_fig.ax_heatmap.set_xlabel("Inference BED loci")
        # Protein names on the LEFT. seaborn puts the row dendrogram there and the tick
        # labels on the right, so swap them: the dendrogram moves to the right of the
        # heatmap (x-inverted so its root still points away from the data) and the labels
        # take the space it vacated.
        _hm = heatmap_fig.ax_heatmap
        _rd = heatmap_fig.ax_row_dendrogram
        _p_hm, _p_rd = _hm.get_position(), _rd.get_position()
        _rd.set_position([_p_hm.x1 + 0.012, _p_rd.y0, _p_rd.width, _p_rd.height])
        _rd.invert_xaxis()
        _hm.yaxis.tick_left()
        _hm.yaxis.set_label_position("left")
        _hm.set_ylabel(f"Proteins  [n] = rank by {rank_desc}")
        plt.setp(_hm.get_yticklabels(), rotation=0, fontsize=8)
        if not args.no_clustering:
            legend_handles = [
                Patch(facecolor=cluster_to_color[cid], edgecolor="none", label=f"C{cid}")
                for cid in cluster_ids_sorted
            ]
            heatmap_fig.ax_heatmap.legend(
                handles=legend_handles,
                title="Locus cluster",
                loc="upper left",
                bbox_to_anchor=(1.0 + (_p_rd.width / max(_p_hm.width, 1e-9)) + 0.05, 1.0),
                frameon=False,
                fontsize=8,
                title_fontsize=9,
            )
        # Default clustermap puts the colourbar top-left, where it now sits on top of the
        # relocated row labels. Park it under the legend on the right instead.
        _p_leg_x = 1.0 + (_p_rd.width / max(_p_hm.width, 1e-9)) + 0.05
        heatmap_fig.ax_cbar.set_position([
            _p_hm.x0 + _p_leg_x * _p_hm.width,
            _p_hm.y0 + 0.05 * _p_hm.height,
            0.015,
            max(0.12, 0.28 * _p_hm.height),
        ])
        heatmap_fig.ax_cbar.tick_params(labelsize=7)
        heatmap_fig.ax_cbar.set_ylabel(cbar_label, fontsize=8)

        heatmap_path = outdir / "binf_support_heatmap.png"
        heatmap_fig.savefig(heatmap_path, dpi=200, bbox_inches="tight")
        plt.close(heatmap_fig.fig)
        print(
            f"Wrote heatmap to: {heatmap_path}"
            + (" (loci ordered by total support; k-means skipped)" if args.no_clustering else " (clustered)")
        )
        # binf_heatmap_clusters.tsv is the ONLY carrier of cluster membership, and both
        # plot_offset_distribution.py --clusters and plot_cluster_metaprofile_at_loci.py
        # --clusters read it. It is therefore written whenever clusters exist, and skipped
        # only when there are none to record.
        all_cluster_labels = None
        if args.no_clustering:
            print("Row clustering skipped (--no-clustering); no binf_heatmap_clusters.tsv written")
        else:
            cluster_sizes = {int(cid): int(np.sum(row_clusters == cid)) for cid in np.unique(row_clusters)}
            print(f"Row cluster sizes (binary k-means, k={k_fit}): {cluster_sizes}")

            # Export per-locus cluster assignment aligned to original inference BED order.
            all_cluster_labels = np.full(n_binf, -1, dtype=np.int64)
            keep_indices = np.flatnonzero(heatmap_keep_mask)
            all_cluster_labels[keep_indices] = row_clusters
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
                    cluster_label = str(int(all_cluster_labels[i])) if passes else "NA"
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
                panel_anchor=args.panel_anchor,
            )

        if args.inspect_protein is not None:
            if inspect_counts_matrix is None:
                available = ", ".join(protein_names)
                raise ValueError(f"Protein '{args.inspect_protein}' not found. Available proteins: {available}")
            nt_offsets = np.arange(-args.window, args.window + 1, dtype=np.int64)
            insp_lo, insp_hi = 0, len(nt_offsets)
            if args.nt_heatmap_window is not None:
                insp_lo = args.window - args.nt_heatmap_window
                insp_hi = args.window + args.nt_heatmap_window + 1
                nt_offsets = nt_offsets[insp_lo:insp_hi]
            inspect_totals = inspect_counts_matrix.sum(axis=1)
            inspect_keep_mask = inspect_totals >= 10
            inspect_counts_filtered = inspect_counts_matrix[inspect_keep_mask, :]
            # Row filter above used the full window on purpose; crop only what is drawn.
            inspect_counts_log = np.log1p(inspect_counts_filtered[:, insp_lo:insp_hi])
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
                # all_cluster_labels is None under --no-clustering; the plotting function
                # already falls back to a single-colour scatter in that case.
                make_tsne_from_binf_summary(
                    cluster_to_color=cluster_to_color,
                    tsv_path=tsv_path,
                    out_path=tsne_path,
                    perplexity=args.tsne_perplexity,
                    random_state=args.tsne_random_state,
                    row_clusters=all_cluster_labels,
                    feature_proteins=cluster_protein_names,
                )


if __name__ == "__main__":
    main()

