# `intersect_inference_bed.py`

Summarize CLIP/iCLIP/eCLIP-style **crosslink (xl) support** around loci from an **inference BED** across *multiple proteins*, producing:

- a **metaprofile plot** (primary output): Gaussian-smoothed **mean** XL file-support signal from \(-window:+window\) (top 10 proteins by default)
- an **optional per-locus summary table** (TSV): per-locus signal shape metrics for each protein
- **per-protein merged XL BEDs** written under `<xldir>/merged/`
- a **clustered heatmap** across all proteins (`binf` rows x proteins columns; values = per-locus total support after logistic scaling)
- a **heatmap cluster assignment table** for downstream cluster-specific metaprofiles
- an optional **single-protein nucleotide heatmap** (`binf` rows x nt columns) via `-i/--inspect-protein`

## Inputs

### Crosslinks directory (`--xldir/-x`)

`--xldir` should point to a directory with **one subdirectory per protein**, e.g.:

```text
xldir/
  PRPF8/
    PRPF8_HepG2_1_R1.genome.xl.bed
    PRPF8_HepG2_2_R1.genome.xl.bed
  FUS/
    ...
```

Each protein directory must contain one or more `*genome.xl.bed` files in **BED6** format:

```text
chrom  start  end  name  score  strand
```

- **score** (column 5) is read but not used directly in merging; merging uses site presence per file
- **strand** (column 6) is required (overlaps are strand-aware)

If `--xldir` itself contains `*genome.xl.bed` (no subdirs), the script treats `--xldir` as a single-protein directory (legacy layout).

**Recommended to use a separate control directory for each protein crosslink profile generated **
example: xl/prpf8 xl/prpf8ctrl xl/fus xl/fusctrl

### Inference BED (`--bed/-b`)

The inference BED must have **at least 6 columns** with strand in column 6:

```text
chrom  start  end  name  score  strand  ...
```

Only the first 6 columns are required; extra columns are ignored.

### Genome sizes (`--genome`)

Used for `bedtools slop` when expanding inference loci by `--window`.

Created by using samtools faidx or cut -f1,2 on reference genome fasta index file (fa.fai file ext.)

## CLI

```bash
python3 scripts/intersect_inference_bed.py \
  -x <xldir> \
  -b <inference.bed> \
  --window 100 \
  --gaussian-sigma 2.0 \
  -i PRPF8 \
  --table \
  -o results/
```

### Options

- **`-x/--xldir`**: xl root directory (required)
- **`-b/--bed`**: inference BED (required)
- **`--window`**: half-window size (default 100). Output offsets run from \(-window..+window\).
- **`--gaussian-sigma`**: sigma parameter for Gaussian metaprofile smoothing (default 2.0)
- **`--cluster-metaprofiles`**: if set, write one metaprofile plot per heatmap cluster (`metaprofile_cluster_C*.png`)
- **`--n-clusters`**: number of k-means clusters on the binarized heatmap row matrix (default 20; capped by number of filtered rows)
- **`--cluster-top-proteins`**: top K XL groups used for heatmap/clustering/tSNE features (default 100)
- **`--metaprofile-top-proteins`**: top K proteins plotted in global/per-cluster metaprofiles (default 15)
- **`-i/--inspect-protein`**: optional protein name for an extra per-nucleotide heatmap for that protein
- **`--skip-merge`**: skip per-protein merge and use direct BED/BED.GZ inputs from `--xldir`
- **`-s/--samplesheet`**: optional TSV (`file`, `group`) used with `--skip-merge`; `file` is resolved relative to `--xldir`
- **`--table`**: if set, write the per-locus summary TSV
- **`--tsne`**: if set, generate tSNE from `binf_summary.tsv` using the same top-K `*_total_overlaps` columns as the heatmap (requires `--table`)
- **`--tsne-perplexity`**: tSNE perplexity (default 30; clipped to valid range)
- **`--tsne-random-state`**: tSNE random seed (default 42)
- **`-o/--outdir`**: directory for plot + TSV (default `results/` in the current working directory)

## What the script does

### 1) Merge xl sites per protein

For each protein directory, all `*genome.xl.bed` files are merged by exact locus and strand.
Each site is tagged by source filename, then grouped so the merged score becomes:

- group key: `(chrom, start, end, strand)`
- aggregation: `count_distinct(file)` (number of xl files containing that exact site+strand)

Output:

```text
<xldir>/merged/<protein>_merged.bed
```

Chromosome names are normalized to `chr*` (e.g. `1` → `chr1`) to match typical BED naming.

### 1b) Alternate input mode (`--skip-merge`)

When `--skip-merge` is used, the script does not merge protein subdirectories.

- with `-s/--samplesheet`: reads a TSV with `file` and `group` columns
  - `file` paths are relative to `--xldir`
  - `group` becomes the protein label in plots/tables
- without `--samplesheet`: uses all BED/BED.GZ files directly in `--xldir`

### 2) Build per-locus signal vectors

Each inference locus is expanded by `--window` (bedtools slop). The script intersects merged xl sites with these windows **strand-aware** (`bedtools intersect -s`).

For each locus, it builds a vector of length \(2*window+1\), where each position stores the **summed file-support score** at that relative offset.

Offsets are **strand-aligned**:

- locus `+`: offset \(=\) `xl_start - locus_start`
- locus `-`: offset \(=\) `-(xl_start - locus_start)` (so + offsets are always in the locus’ 5'→3' direction)

## Outputs

### Metaprofile plot (always)

File: `<outdir>/metaprofile.png`

For each protein:

- compute `total_overlaps = sum(vector)` for each locus
- compute the **mean** support vector across all loci (average by number of input `binf` regions)
- smooth with a **Gaussian kernel** controlled by `--gaussian-sigma`
- rank proteins by total smoothed metaprofile signal and plot only the **top K** from `--metaprofile-top-proteins` (default 15)
- place legend on the right side of the figure

### Clustered heatmap (always)

File: `<outdir>/binf_support_heatmap.png`

- rows: `binf` loci
- columns: the **top K** XL groups by global total signal from `--cluster-top-proteins` (default 100), not the full protein list
- values: per-locus total support (`total_overlaps`) transformed by logistic scaling
- pre-filter rows: keep loci with `sum(total_overlaps across *all* XL groups) >= 40` (filter uses the full table; heatmap columns are still top-K only)
- **row clustering**: build a **binary** matrix over the top-K proteins (`1` if support `> 0` at that locus/protein, else `0`), then **k-means** (`--n-clusters`, default 20) on those binary rows
- **row order** (no row dendrogram): sort by cluster id (ascending), then by total crosslink support across all proteins (descending) within each cluster
- **column clustering only**: cosine distance + average linkage between protein columns (computed on the scaled matrix before row reorder; row order does not change column vectors for linkage)
- heatmap **values** shown are still **logistic-scaled** continuous totals (viridis)

### Heatmap cluster assignments (always)

File: `<outdir>/binf_heatmap_clusters.tsv`

- one row per input `binf` locus (same order as summary table)
- columns:
  - `binf_chr_start_end`
  - `chrom`, `start`, `end`
  - `row_sum_support` (sum of `total_overlaps` across **all** XL groups)
  - `passes_heatmap_filter` (`True`/`False`, threshold `>=40` on sum across all proteins)
  - `heatmap_cluster`: `NA` if row fails heatmap filter; otherwise integer k-means cluster id (`1..k`)

### Cluster metaprofiles (optional)

Enabled by `--cluster-metaprofiles`.

- one metaprofile plot per heatmap row cluster
- files: `<outdir>/metaprofile_cluster_C1.png`, `<outdir>/metaprofile_cluster_C2.png`, ...
- for each cluster: mean support profile is computed only from loci assigned to that cluster

### tSNE from summary table (optional)

Enabled by `--table --tsne`.

File: `<outdir>/binf_summary_tsne.png`

- input features: `*_total_overlaps` columns for the **same top-K proteins** as the heatmap (`--cluster-top-proteins`)
- feature transform: same logistic scaling used for the global heatmap
- one point per `binf` row
- points: colored by k-means cluster id (`1..k`); rows not in heatmap filter are light gray
- requires `scikit-learn` in the environment

### Single-protein nucleotide heatmap (optional)

Enabled by `-i/--inspect-protein`.

File: `<outdir>/binf_<protein>_nt_support_heatmap.png`

- rows: `binf` loci
- columns: nucleotide positions from `-window..+window`
- values: merged support score at each relative nucleotide offset for the selected protein
- pre-filter rows: keep loci with `sum(across nt positions) >= 10` for the selected protein
- hierarchical clustering on rows only; nucleotide columns stay in genomic order (`-window..+window`)

### Summary table (optional)

Enabled by `--table`.

File: `<outdir>/binf_summary.tsv`

- **Column 1**: `binf_chr_start_end` formatted as `chr_start_end`
- For each protein, five columns are added:
  - `<protein>_total_overlaps`
  - `<protein>_variance`
  - `<protein>_pearson_median_skew`
  - `<protein>_kurtosis_excess`
  - `<protein>_max_binding_offset`

These metrics are computed from the per-locus vector across \(-window..+window\):

- **total_overlaps**: \(\sum\) of xl scores across the full vector
- **variance**: variance of xl scores across the full vector
- **pearson_median_skew**: Pearson’s median skewness

  \[
  3 \cdot \frac{\text{mean} - \text{median}}{\text{std}}
  \]

  (defined as 0 when `std == 0`).

- **kurtosis_excess**: Fisher excess kurtosis; Positive values indicate leptokurtic distribution with strong tailedness while negative values are strongly platykurtic

  \[
  \frac{\mu_4}{\sigma^4} - 3
  \]

- **max_binding_offset**: offset of maximal binding based on **sliding 5-nt window sums**
  - compute all 5-nt window sums along the vector
  - take the window with the maximum sum
  - report the **center** position of that window as an offset

## Notes

- Input coordinates can repeat (duplicate `chr/start/end`); the script keeps **one row per input BED line** and handles duplicates during intersection.
- Runtime is dominated by the `bedtools groupby` merge and `bedtools intersect` steps for large XL datasets.
- `-i/--inspect-protein` must match a protein directory name under `--xldir`.

