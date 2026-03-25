# RBPeek
merge crosslinks files and intersect with BED for broad binding analysis and metaprofiling
# `intersect_decoy.py`

Summarize CLIP/iCLIP/eCLIP-style **crosslink (XL) signal** around **decoy splice-site loci** across *multiple proteins*, producing:

- a **metaprofile plot** (primary output): smoothed aggregate XL signal from \(-window..+window\) around the decoy site
- an **optional per-decoy summary table** (TSV): per-decoy overlap/signal shape metrics for each protein
- **per-protein merged XL BEDs** for reuse and reproducibility

## Inputs

### Crosslinks directory (`--xldir/-x`)

`--xldir` should point to a directory with **one subdirectory per protein**, e.g.:

```text
xldir/
  PRPF8/
    PRPF8_HepG2_1_R1.genome.xl.bed
    PRPF8_HepG2_2_R1.genome.xl.bed
  HNRNPH1/
    ...
```

Each protein directory must contain one or more `*genome.xl.bed` files in **BED6** format:

```text
chrom  start  end  name  score  strand
```

- **score** (column 5) must be numeric (XL count/weight per site).
- **strand** (column 6) is used and overlaps are computed **strand-aware**.

If `--xldir` itself contains `*genome.xl.bed` (no subdirs), the script treats `--xldir` as a single-protein directory (legacy layout).

### Decoys BED (`--decoys/-d`)

Decoy BED must have **at least 6 columns** with strand in column 6:

```text
chrom  start  end  name  score  strand  ...
```

Only the first 6 columns are required; extra columns are ignored.

### Genome sizes (`--genome`)

Used for `bedtools slop` when expanding decoys by `--window`.

Default: `advbfx/reference/genomes/Gencode49/genome.sizes`

## CLI

```bash
python3 intronretention/hnRNPH1_IR_MAF/scripts/intersect_decoy.py \
  -x <xldir> \
  -d <decoys.bed> \
  --window 50 \
  --min-sum 25 \
  --table \
  -o results/
```

### Options

- **`-x/--xldir`**: XL root directory (required)
- **`-d/--decoys`**: decoy BED (required)
- **`--window`**: half-window size (default 50). Output offsets run from \(-window..+window\).
- **`--min-sum`**: minimum summed XL signal used to *filter decoys into the metaprofile only* (default 5)
- **`--table`**: if set, write the per-decoy summary TSV
- **`-o/--outdir`**: directory for plot + TSV (default `results/` in the current working directory)

## What the script does

### 1) Merge XL sites per protein

For each protein directory, all `*genome.xl.bed` files are merged by exact locus and strand:

- group key: `(chrom, start, end, strand)`
- aggregation: `sum(score)` (score is BED column 5)

The merged file is written to:

```text
<xldir>/merged/<protein>_merged.bed
```

Chromosome names are normalized to `chr*` (e.g. `1` → `chr1`) to match decoy inputs.

### 2) Build per-decoy signal vectors

Each decoy is expanded by `--window` (bedtools slop). The script intersects merged XL sites with these windows **strand-aware** (`bedtools intersect -s`).

For each decoy, it builds a vector of length \(2*window+1\), where each position stores the **summed XL score** at that relative offset.

Offsets are **decoy-strand aligned**:

- decoy `+`: offset \(=\) `xl_start - decoy_start`
- decoy `-`: offset \(=\) `-(xl_start - decoy_start)` (so + offsets are always in the decoy’s 5'→3' direction)

## Outputs

### Metaprofile plot (always)

File: `<outdir>/metaprofile.png`

For each protein:

- compute `total_overlaps = sum(vector)` for each decoy
- keep only decoys where `total_overlaps >= --min-sum`
- sum vectors across kept decoys to get an aggregate profile
- smooth with a **centered rolling sum over 5 nt**
- plot all proteins as separate curves on the same axes

### Summary table (optional)

Enabled by `--table`.

File: `<outdir>/decoy_summary.tsv`

- **Column 1**: `decoy_chr_start_end` formatted as `chr_start_end`
- For each protein, five columns are added:
  - `<protein>_total_overlaps`
  - `<protein>_variance`
  - `<protein>_pearson_median_skew`
  - `<protein>_kurtosis_excess`
  - `<protein>_max_binding_offset`

These metrics are computed from the per-decoy vector across \(-window..+window\):

- **total_overlaps**: \(\sum\) of XL scores across the full vector.
- **variance**: variance of XL scores across the full vector (\(-window..+window\)).
- **pearson_median_skew**: Pearson’s *median* skewness

  \[
  3 \cdot \frac{\text{mean} - \text{median}}{\text{std}}
  \]

  (defined as 0 when `std == 0`).

- **kurtosis_excess**: Fisher **excess kurtosis**

  \[
  \frac{\mu_4}{\sigma^4} - 3
  \]

  Normal-like tails correspond to ~0; larger values indicate a more “spiky/outlier-heavy” profile.

- **max_binding_offset**: offset of maximal binding based on **sliding 5-nt window sums**:
  - compute all 5-nt window sums along the vector
  - take the window with the maximum sum
  - report the **center** position of that window as an offset

## Notes / gotchas

- Decoy coordinates can repeat (duplicate `chr/start/end`); the script keeps **one row per input decoy line** in the TSV and handles duplicates during intersection.
- Runtime is dominated by the `bedtools groupby` merge and `bedtools intersect` steps for large XL datasets.

