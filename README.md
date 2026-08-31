# `intersect_inference_bed.py`

Summarize CLIP/iCLIP/eCLIP-style **crosslink (xl) support** around loci from an **inference BED** across *multiple proteins*, producing:

- a **metaprofile plot** (primary output): Gaussian-smoothed **mean** peak support across all loci from \(-window:+window\), for the top 10 proteins of the **same ranking that picks the heatmap's columns**, with a right-hand axis giving the equivalent total summed score
- an **optional per-locus summary table** (TSV): per-locus signal shape metrics for each protein
- **per-protein merged XL BEDs** written under `<xldir>/merged/` (merge mode only; not under `--skip-merge`)
- a **clustered heatmap** across the top-K proteins (`binf` rows x proteins columns; values = per-locus total support, scaled per `--heatmap-scale`)
- a **ranked panel-sample table** (`protein_ranking.tsv`): every column ranked by mean peak support, carrying its central-enrichment score alongside
- a **heatmap cluster assignment table** for downstream cluster-specific metaprofiles (unless `--no-clustering`)
- an optional **single-protein nucleotide heatmap** (`binf` rows x nt columns) via `-i/--inspect-protein`

This repo also carries the analysis wrappers built around the tool. See
[Repository layout](#repository-layout) for what each script does, and
[`THRAP3/README.md`](THRAP3/README.md) and [`Centrosome/README.md`](Centrosome/README.md)
for the two analyses and their findings.

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
- **`--panel-anchor`**: `start` (default) or `midpoint` — which point of each `--xldir` interval carries its score. Use `midpoint` whenever the panel holds **peaks** rather than 1 nt crosslink sites; see [Panel anchor](#panel-anchor-start-vs-midpoint) below.
- **`--nt-heatmap-window`**: crop the `-i/--inspect-protein` nucleotide heatmap to +/- this many nt for **display only** (default: full `--window`). `--window` still governs what is counted, so all statistics, the metaprofile, the enrichment ranking and the summary table are unchanged.
- **`--heatmap-min-support`**: minimum summed support across **all** xl groups for a locus to enter the heatmap/clustering (default 50). Scales with panel width, so raise it for wide panels — see the note under the clustered heatmap.
- **`--heatmap-min-support-percentile`**: express that filter as "drop the bottom P% of loci" instead, and always drop zero-support loci. Overrides `--heatmap-min-support`. Use it whenever two runs will be compared — see [the row filter](#the-heatmap-row-filter).
- **`--heatmap-scale`**: `logistic` (default) or `percentile` — colour scaling for the clustered heatmap; see [Heatmap colour scaling](#heatmap-colour-scaling) below.
- **`--heatmap-scale-percentile`**: percentile of non-zero values mapped to the top of the colour range under `--heatmap-scale percentile` (default 99.0)
- **`--protein-select`**: `total` (default), `enrichment` (recommended) or `centrality` — how the top-K columns are chosen; see [Protein selection](#protein-selection-total-enrichment-centrality) below.
- **`--enrichment-window`**: half-width (nt) counted as "on the locus" by `enrichment` ranking (default 5)
- **`--protein-min-loci`**: minimum loci with signal for a protein to be eligible under `enrichment`/`centrality` (default 500)
- **`--protein-min-loci-frac`**: express that gate as a **fraction of the inference loci** instead, so it scales with the set. Overrides `--protein-min-loci`. Use it whenever two runs will be compared — see [the eligibility gate](#the-eligibility-gate).
- **`--centrality-sigma`**: width (nt) of the Gaussian template for centrality ranking (default 5.0)
- **`--centrality-min-total`**: minimum summed `total_overlaps` for a protein to be eligible under centrality ranking (default 100)
- **`--exclude-groups`**: comma-separated substrings; panel columns whose group name contains one are dropped before anything is computed (e.g. `--exclude-groups PARCLIP`)
- **`--no-clustering`**: skip k-means row clustering entirely. The heatmap is still written, with loci ordered by total support; the tSNE becomes a single-colour scatter; `binf_heatmap_clusters.tsv` is not written. Incompatible with `--cluster-metaprofiles`.
- **`--cluster-metaprofiles`**: if set, write one metaprofile plot per heatmap cluster (`metaprofile_cluster_C*.png`)
- **`--n-clusters`**: number of k-means clusters on the binarized heatmap row matrix (default 20; capped by number of filtered rows)
- **`--cluster-top-proteins`**: top K XL groups used for heatmap/clustering/tSNE features (default 60)
- **`--metaprofile-top-proteins`**: top K proteins plotted in global/per-cluster metaprofiles (default 10). These are the first K of the **same** `--protein-select` ranking that picks the heatmap's columns, so the metaprofile is always a subset of the heatmap; capped at `--cluster-top-proteins`.
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

### 1c) Chromosome naming is harmonised automatically

Merge mode normalises panel chromosome names to `chr*`; `--skip-merge` does not, so a panel
file named Ensembl-style (`1`) against a `chr1` inference BED used to produce **zero**
overlaps. That failure is silent — the run completes and the column simply reads as an RBP
that binds nothing. Five columns were lost that way before this was fixed.

The script now compares each panel file's first data chromosome against the inference BED's
convention and rewrites only the mismatched ones into its temp directory, logging:

```text
Chromosome naming: rewrote 5 panel column(s) to match the inference BED (chr-prefixed style): ...
```

A correctly-named panel costs one line read per file and is otherwise untouched. Panel files
with no data rows are warned about, since they look identical in the output.

### 2) Build per-locus signal vectors

Each inference locus is expanded by `--window` (bedtools slop). The script intersects merged xl sites with these windows **strand-aware** (`bedtools intersect -s`).

For each locus, it builds a vector of length \(2*window+1\), where each position stores the **summed file-support score** at that relative offset.

Offsets are **strand-aligned**:

- locus `+`: offset \(=\) `xl_start - locus_start`
- locus `-`: offset \(=\) `-(xl_start - locus_start)` (so + offsets are always in the locus’ 5'→3' direction)

### Panel anchor: `start` vs `midpoint`

A panel interval's **entire score lands on one offset**, it is not spread across its width.
Which offset that is comes from `--panel-anchor`.

For the 1 nt crosslink sites this script was designed around, the two settings are
identical (`start == midpoint` when `end == start + 1`). They diverge as soon as the panel
holds **peaks**, and the divergence is not a harmless shift:

- a peak of width `w` centred on a `+` strand locus scores at **`-w/2`**
- the strand flip sends the same peak on a `-` strand locus to **`+w/2`**, because a peak's
  genomic start is its 3' end there

So real central binding is split into a spurious **`+/-w/2` doublet with a hole at 0**.
Measured on a panel of Clippy peaks (mean width ~11 bp) against 29,018 centred loci:

| `--panel-anchor` | `+` strand median | `-` strand median | separation |
|---|---:|---:|---:|
| `start` | -5.0 | +5.0 | **10.0 nt** |
| `midpoint` | 0.0 | 0.0 | **0.0 nt** |

**Use `midpoint` for any peak-based panel.** It is a no-op for crosslink input, so it is
safe to leave on. `start` remains the default only so existing runs reproduce byte for byte.

Note that `max_binding_offset` in the summary table centres a 5 nt sliding window, so it
carries its own quantisation of a couple of nt independently of this setting.

### Heatmap colour scaling

`logistic` (default) centres on the matrix **median**. The locus x protein matrix is sparse,
so that median is ~0 and **every empty cell maps to exactly 0.5** — mid-palette. The whole
range below 0.5 goes unused and the colourbar starts at 0.5, which is why a sparse run looks
uniformly flat no matter how the row filter is set.

`percentile` applies `log1p`, then scales against the given percentile of the **non-zero**
values and clips. Empty cells stay at 0, so the full palette carries signal. On a sparse test
matrix (65% of cells empty):

| scaling | range | empty cells render at | median non-zero cell |
|---|---|---:|---:|
| `logistic` | 0.50 – 1.00 | 0.50 | 0.67 |
| `percentile` (log1p, 99th) | 0.00 – 1.00 | 0.00 | 0.53 |

The `log1p` step is what makes this usable rather than merely correct. Support counts are
heavy-tailed — non-zero median 10 against a maximum of 333 on that matrix — so scaling raw
values against the 99th percentile put the median cell at **0.11**, *darker* than the
logistic scaling it replaces. After `log1p` the same settings put it at 0.53.

Two things this also affects:

- the **column dendrogram**, built from cosine distances on the scaled matrix. Under
  `logistic` every column carries a large constant 0.5 component from its empty cells, which
  compresses the distances between them; under `percentile` the distances reflect actual
  co-occurrence.
- **not** the row clusters, which come from `(matrix > 0)` binarisation and are unaffected by
  any colour scaling.

Raising `--heatmap-min-support` is *not* an alternative fix. It culls rows rather than
rescaling colour, and because row support correlates with inference-BED reproducibility it
biases which loci survive.

### The heatmap row filter

Loci enter the heatmap and clustering only if their support clears a threshold. Support here
is the **summed score column of every overlapping panel interval, across every column and the
whole window** — so it has no intrinsic scale. It rises with panel width, with sequencing
depth, and with how bound the region is, and none of those are properties of the locus.

That is why an absolute `--heatmap-min-support` cannot be shared by two runs. Measured on the
THRAP3 exonic and intronic subsets, median row support was **1,594** and **296** — a 5.4x
gap — so a flat threshold of 200 dropped **11.3%** of exonic loci and **43.5%** of intronic
ones. The intronic heatmap was a top-56% subset of its loci while the exonic one was a
top-89% subset, and the two were being read side by side.

`--heatmap-min-support-percentile P` fixes that by dropping the bottom P% instead, and is the
right choice whenever two runs will be compared. It **always** drops zero-support loci as
well, whatever the percentile resolves to. That clause is load-bearing: 18.6% of the THRAP3
intronic loci have no support from any panel column, so the 10th percentile of that set is
literally 0 and a bare percentile would keep ~750 blank heatmap rows. At `P=10` the exonic
set resolves to a threshold of 166 and loses 10.0%; the intronic set resolves to 0, the zero
clause fires, and it loses 18.6%. The residual asymmetry is then real emptiness rather than
an artefact of scale.

The same argument applies to `--protein-min-loci` vs `--protein-min-loci-frac` for the
column-eligibility gate — see [the eligibility gate](#the-eligibility-gate).

### Protein selection: `total`, `enrichment`, `centrality`

`--cluster-top-proteins K` picks which panel columns reach the heatmap, clustering, the tSNE
and — as its first `--metaprofile-top-proteins` entries — the metaprofile. `--protein-select`
decides *how* they are picked.

`total` (default) ranks by summed `total_overlaps`. That is a **depth-biased** measure: a
deeply sequenced, peak-rich dataset scores highly whether or not its binding has anything
to do with the inference loci.

`centrality` ranks by Pearson r between each protein's mean profile and a Gaussian of width
`--centrality-sigma` centred on the locus. Pearson r is **scale-free**, so a shallow dataset
with sharply centred binding can outrank a deep one with a flat profile.

`enrichment` (**recommended**) is also depth-free but works **per locus**: for each protein,
the fraction of its signal-bearing loci whose `max_binding_offset` falls within
`--enrichment-window`. Prefer it over `centrality` — see below for why.

#### Why `enrichment` over `centrality`

Once the panel is centred with `--panel-anchor midpoint`, nearly every mean profile peaks at
0, so correlating against a centred template stops discriminating on **position** and starts
discriminating on how tidy a dataset's **baseline** is. That systematically favours assays
with low background. Measured on a 297-column run:

| | `centrality` | `enrichment` |
|---|---|---|
| assay mix of top 20 | 10 PAR-CLIP / 9 eCLIP / 1 iCLIP | 19 eCLIP / 1 iCLIP |
| vs panel composition (75% eCLIP, 23% PAR-CLIP) | PAR-CLIP **2.2x** enriched | matches |
| BCLAF1 (known complex partner of the bait) | **absent from top 20** | 16th |
| top-ranked columns | ORF1/L1RE1-PARCLIP (LINE-1, repeat-derived) | both HNRNPC datasets |

The two rankings shared **zero** of their top 20 on the same data.

`--protein-min-loci` is load-bearing for `enrichment`: the fraction is trivially 1.0 for a
protein whose signal touches a single locus. Ungated, the top of that ranking was a PAR-CLIP
column with one locus and a summed score of 5.

Validated against two decoys built from a real replicate at matched depth — one displaced
40 nt, one with positions jittered +/-80 nt:

| column | summed overlaps | rank by `total` | centrality r | rank by `centrality` |
|---|---:|---:|---:|---:|
| genuine replicate A | 160,016 | 2 | +0.813 | **1** |
| genuine replicate B | 164,279 | 1 | +0.809 | **2** |
| genuine replicate C | 104,214 | 5 | +0.808 | **3** |
| genuine replicate D | 43,386 | 6 | +0.796 | **4** |
| decoy, jittered | 157,163 | 4 | +0.292 | 5 |
| decoy, shifted 40 nt | 158,592 | 3 | **-0.003** | 6 |

Under `total` both decoys outrank two genuine replicates on depth alone. Under `centrality`
every genuine replicate outranks both decoys, including one at a quarter of their depth.
The jittered decoy retains r=+0.292 because jitter is itself centred, leaving broad central
enrichment inside the window — a broad hump is genuinely less centred, not a scoring flaw.

Proteins below `--centrality-min-total` are excluded before ranking, since a profile built
from a handful of overlaps can score a near perfect correlation off a single spike at 0.
Flat profiles are excluded too, Pearson r being undefined at zero variance. Selected
columns are logged with both their r and their total, so a high-r/low-signal column is
visible rather than silently shaping the heatmap.

#### The eligibility gate

`enrichment` scores a **fraction with a variable denominator** — of the loci where a protein
has any signal, how many peak within `--enrichment-window`. A protein with signal at one
locus that happens to peak on target scores a perfect **1.000**, beating a protein at 0.380
built on 7,968 loci. Ungated, the top of one real ranking was a PAR-CLIP column with **one**
locus and a summed score of 5.

Two gates apply, and a protein must clear both or it is dropped from the ranking entirely:

- `--centrality-min-total` (default 100) — enough total signal
- `--protein-min-loci` (default 500) — enough loci for the fraction to mean anything

Precision of the fraction, near 0.30:

| loci | standard error |
|---:|---|
| 100 | ±0.046 |
| 300 | ±0.026 |
| 500 | ±0.021 |
| 5,000 | ±0.006 |

**`--protein-min-loci` is an absolute count, so it does not scale with the inference set.**
Pass **`--protein-min-loci-frac`** instead whenever two runs will be compared — it expresses
the gate as a fraction of that run's loci, so both are gated identically in relative terms,
and the run log prints the value it resolved to.

Comparing two runs of different size therefore gates them unequally: 500 demanded signal at
1.7% of loci on a 29,018-locus set but 5.8% on an 8,666-locus one, excluding 38 versus 90 of
302 proteins. A protein can then be ranked in one run and invisible in the other purely from
set size — and **exclusion is indistinguishable from absence** in the output. When comparing
runs, size-match the gate (`500 x n_small/n_large`) and check the `N of M proteins excluded`
line in both logs.

## Outputs

### Metaprofile plot (always)

File: `<outdir>/metaprofile.png`

For each protein:

- compute `total_overlaps = sum(vector)` for each locus
- compute the **mean** support vector across all loci (`counts.mean(axis=0)`) — note the
  denominator is **every** locus, including those where the protein has no signal at all, so
  a curve is diluted by non-binding loci rather than describing the sites it does bind
- smooth with a **Gaussian kernel** controlled by `--gaussian-sigma`
- plot the **top K** from `--metaprofile-top-proteins` (default 10), taken from the **same
  `--protein-select` ranking that picks the heatmap's columns**. The metaprofile is therefore
  always a subset of the heatmap, and the two figures never disagree about which proteins
  matter. (Before this, the metaprofile ranked by summed profile signal — pure depth — while
  the heatmap ranked by enrichment, so they routinely showed different proteins.)
- **left axis**: mean peak support across all loci. **Right axis**: the same curve multiplied
  by the locus count, i.e. total summed panel score. Because the locus count is one constant,
  the two axes are the same curve at two scales and agree pixel for pixel — the only
  normalisation that can honestly share an axis. Anything per-protein (dividing each curve by
  its own maximum, say) reorders the curves and needs its own panel.
- legend on the right, each entry annotated with that protein's grand total

### Clustered heatmap (always)

File: `<outdir>/binf_support_heatmap.png`

- rows: `binf` loci
- columns: the **top K** XL groups from `--cluster-top-proteins` (default 60), ranked per `--protein-select`, not the full protein list
- values: per-locus total support (`total_overlaps`) transformed by logistic scaling
- pre-filter rows: see [the heatmap row filter](#the-heatmap-row-filter) below
- **row clustering**: build a **binary** matrix over the top-K proteins (`1` if support `> 0` at that locus/protein, else `0`), then **k-means** (`--n-clusters`, default 20) on those binary rows. Skipped entirely under `--no-clustering`.
- **row order** (no row dendrogram): sort by cluster id (ascending), then by total support across all proteins (descending) within each cluster. Under `--no-clustering`, by total support alone.
- **column clustering only**: cosine distance + average linkage between protein columns (computed on the scaled matrix before row reorder; row order does not change column vectors for linkage)
- heatmap **values** shown are still **logistic-scaled** continuous totals (viridis)

### Ranked panel samples (always)

File: `<outdir>/protein_ranking.tsv`

One row per panel column, sorted by `mean_peak_support` descending. Columns:

| column | meaning |
|---|---|
| `rank` | position by `mean_peak_support` |
| `sample` | samplesheet `group` |
| `mean_peak_support` | `total_peak_support / n_loci` — the metaprofile's left axis, integrated |
| `total_peak_support` | summed panel score across every locus and offset |
| `loci_with_signal`, `frac_loci_with_signal` | how much of the locus set this column touches at all |
| `<select>_score` | the `--protein-select` score — `enrichment_score` is the fraction of signal-bearing loci peaking within `--enrichment-window`; `NA` if the column was gated out |
| `<select>_rank` | rank by that score; `NA` if gated out |
| `eligible` | passed `--protein-min-loci*` and `--centrality-min-total` |
| `selected_for_heatmap` | reached the heatmap's top K |

This table carries **both** competing answers to "which RBP co-occupies these sites" side by
side: depth (`mean_peak_support`) and central enrichment (`frac centred`). A column that is
deep but binds 30 nt away ranks 1st on the former and last on the latter, and the table makes
that visible instead of forcing the choice up front.

### Heatmap cluster assignments (unless `--no-clustering`)

File: `<outdir>/binf_heatmap_clusters.tsv`

- one row per input `binf` locus (same order as summary table)
- columns:
  - `binf_chr_start_end`
  - `chrom`, `start`, `end`
  - `row_sum_support` (sum of `total_overlaps` across **all** XL groups)
  - `passes_heatmap_filter` (`True`/`False`, per `--heatmap-min-support` / `--heatmap-min-support-percentile`)
  - `heatmap_cluster`: `NA` if row fails heatmap filter; otherwise integer k-means cluster id (`1..k`)

This file is the **only** carrier of cluster membership, and both
`scripts/plot_offset_distribution.py --clusters` and
`scripts/plot_cluster_metaprofile_at_loci.py --clusters` read it. `--no-clustering` suppresses
it, so those two scripts need a run that kept clustering.

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
- points: colored by k-means cluster id (`1..k`), in the **same** hues the heatmap's cluster bar uses, so `C3` here is `C3` there; rows not in the heatmap filter are light gray. Under `--no-clustering` this becomes a single-colour scatter.
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

## Repository layout

`intersect_inference_bed.py` is the analysis engine. Everything else prepares its inputs or
interprets its outputs.

### Building an inference BED

| script | purpose |
|---|---|
| `build_inference_bed_from_peaks.py` | **General.** Replicate Clippy peak calls -> inference BED. Normalises chromosome names, merges strand-aware, keeps regions supported by `--min-reps` replicates, collapses each to a 1 nt midpoint anchor. `--subtract` drops regions overlapping control peaks, for proximity-labelling baits. |
| `build_thrap3_inference_bed.py` | The THRAP3-specific original, kept for reproducibility of that analysis. Superseded by the general script above. |
| `merge_replicate_peaks.py` | Pools replicate peak BEDs by SUMMING scores at identical chrom/start/end/strand, then reformats to BED6. Exact-match grouping, so it pools signal rather than assessing reproducibility - for the latter use the overlap-based merge above. |
| `split_inference_bed_by_region.py` | Splits an inference BED into **exonic** and **intronic** subsets from a GTF. Strand-aware, exon-priority (all transcripts' exons merged first), so the two sets are disjoint by construction. Defaults to GENCODE v39. |

Why each transformation is needed is documented in the scripts themselves; the short version
is that all three failures they prevent are **silent** — an unnormalised BED intersects
nothing, an unfiltered union follows sequencing depth, and a locus wider than 1 nt smears the
metaprofile by its own width.

### Panel management

| script | purpose |
|---|---|
| `add_new_peaks_to_samplesheet.py` | Discovers `*_Peaks.bed` and appends rows to `RBPeekSamplesheet.tsv`, with paths relative to `--xldir`. Idempotent, supports `--dry-run`. Skips patterns in `--skip` (default `TRA2A,_Mm,B_cells,Bcells`) — see the script for why each is there. |
| `run_clippy_new_samples.sbatch` | Calls Clippy peaks for samples that only exist as crosslink BEDs, with the parameters the panel was built with. |

### Interpreting results

| script | purpose |
|---|---|
| `compare_binding_frequency.py` | Compares how often each panel protein binds two different locus sets, with a two-proportion z-test. Use when positional metrics are uninformative — notably for proximity labelling, where peak position reflects distance from the bait rather than an RBP footprint. |
| `plot_cluster_metaprofile_at_loci.py` | True metaprofile for one heatmap cluster's proteins, computed around a DIFFERENT locus set. Recomputes the counts for just those proteins - the per-nucleotide vector is not in `binf_summary.tsv`, which keeps only five per-locus summary statistics. Imports the offset arithmetic from `intersect_inference_bed.py` so slop, strand handling and `--panel-anchor` match a full run. |
| `plot_offset_distribution.py` | Lighter alternative when a rerun is not wanted: distribution of each locus's single strongest offset, straight from `binf_summary.tsv`. Not a metaprofile. |
| `build_gene_matrix_from_summaries.py` | Flow `*.summary_gene.tsv` files -> gene x sample counts and RPKM matrices. Gene length is in the input, so no GTF is needed; the library denominator is every row's cDNA including intergenic, which is what the file's own `cDNA %` is computed against. |
| `plot_expression_heatmap.py` | Clustered heatmap of normalised expression per gene across samples, for a gene list (`.xlsx` first column, or plain text). Computes RPKM from raw counts with `--gtf`, or takes `--already-normalised`. Reports unmatched genes rather than dropping them silently, including Excel's date-mangled `MARCH`/`SEPT` symbols. |

### Analysis runners

| sbatch | what it runs |
|---|---|
| `run_intersect_inference_bed.sbatch` | Original decoys/splice-site run. |
| `run_thrap3_intersect.sbatch` | THRAP3 against the full panel. |
| `run_thrap3_region.sbatch <exonic\|intronic>` | THRAP3 split by transcript region; identical settings so divergence is region, not parameters. |
| `run_centro_intersect.sbatch` | Centrosome apex CLIP, uncontrolled. |
| `run_centro_controlled.sbatch <specific\|ntcontrol>` | Centrosome with APEX controls applied. Run **both**; neither is interpretable alone. |

### Python version

`intersect_inference_bed.py` needs the `rbpeek` conda env (Python 3.12). The helper scripts
above are written to run under **Python 3.6** as well, since they are typically invoked by
hand on an HPC login node where `python3` is the system interpreter.

## Notes

- Input coordinates can repeat (duplicate `chr/start/end`); the script keeps **one row per input BED line** and handles duplicates during intersection.
- Runtime is dominated by the `bedtools groupby` merge and `bedtools intersect` steps for large XL datasets.
- `-i/--inspect-protein` matches a **protein directory name** under `--xldir` in merge mode, or the **`group` label** from the samplesheet under `--skip-merge`. The error message on a mismatch lists the available names.

