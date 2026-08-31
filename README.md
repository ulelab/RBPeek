# `intersect_inference_bed.py`

Summarise CLIP/iCLIP/eCLIP peak support around the loci of an **inference BED** across a
panel of samples, answering "which RBPs co-occupy these sites, and with what profile shape".

Every run writes five things:

| file | what it is |
|---|---|
| `metaprofile.png` | Gaussian-smoothed mean peak support vs offset, for the top `--top-proteins` samples, with a right-hand axis giving the equivalent total summed score |
| `binf_support_heatmap.png` | loci x samples, values = per-locus total support |
| `protein_ranking.tsv` | every panel sample with **both** rankings — depth and central enrichment — side by side |
| `binf_summary.tsv` | per-locus shape statistics for every sample |
| `binf_summary_tsne.png` | with `--tsne` |

Passing `-n/--n-clusters` adds k-means groups over the loci, `binf_heatmap_clusters.tsv`, and
one metaprofile per cluster.

This repo also carries the analysis wrappers built around the tool. See
[Repository layout](#repository-layout) for what each script does, and
[`THRAP3/README.md`](THRAP3/README.md) and [`Centrosome/README.md`](Centrosome/README.md)
for the two analyses and their findings.

## Inputs

### Samplesheet (`-s/--samplesheet`) and panel root (`-x/--xldir`)

A TSV with two columns:

```text
file	group
eCLIP-Clippy/HepG2-BCLAF1-merged.xl_..._Peaks.bed	HepG2-BCLAF1-eCLIP
```

`file` is resolved relative to `--xldir`; `group` becomes the sample label in every plot and
table. Files may be BED6+ or `.bed.gz`. A missing file is a hard error, not a skipped row.

### Inference BED (`-b/--bed`)

BED6+, strand in column 6; extra columns are ignored. Loci should be **1 nt** — the script
deposits a panel interval's whole score at one offset, so a wider locus smears the
metaprofile by its own width. `build_inference_bed_from_peaks.py` produces 1 nt anchors.

Coordinates may repeat; the script keeps one row per input line and handles duplicates
during intersection.

### Genome sizes (`--genome`)

Used for `bedtools slop` when expanding loci by `--window`. Make it with `samtools faidx`
then `cut -f1,2` on the `.fa.fai`.

## CLI

```bash
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b THRAP3/THRAP3_merged_min2rep_anchors.bed \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  -o results/thrap3 \
  --panel-anchor midpoint \
  --tsne
```

### Options

**Required**

- **`-x/--xldir`** — root the samplesheet's `file` paths resolve against
- **`-b/--bed`** — inference BED
- **`-s/--samplesheet`** — TSV with `file` and `group` columns

**Counting**

- **`--window`** — half-window in bp (default 100); offsets run `-window..+window`
- **`--panel-anchor`** — `start` (default) or `midpoint`. **Use `midpoint` for any panel of peaks**; see [Panel anchor](#panel-anchor-start-vs-midpoint).
- **`--gaussian-sigma`** — metaprofile smoothing sigma (default 2.0)
- **`--genome`** — genome sizes file

**Choosing what to plot**

- **`--protein-select`** — `enrichment` (default) or `total`; see [Protein selection](#protein-selection-total-vs-enrichment)
- **`--enrichment-window`** — half-width (nt) counted as "on the locus" (default 5)
- **`--top-proteins`** — how many samples reach the heatmap, metaprofile and tSNE (default 20). One number for all three, so the figures always show the same set.

**Heatmap**

- **`--support-pct`** — drop the bottom P% of loci by summed support (default 10); zero-support loci are always dropped. See [the row filter](#the-heatmap-row-filter).
- **`--heatmap-scale`** — `percentile` (default) or `logistic`; see [colour scaling](#heatmap-colour-scaling)
- **`--heatmap-scale-percentile`** — non-zero percentile mapped to the top of the range (default 99)

**Optional extras**

- **`-n/--n-clusters`** — k-means the loci into this many groups. Omitted by default: loci are ordered by total support. Passing it also writes `binf_heatmap_clusters.tsv` and one metaprofile per cluster, and colours the tSNE.
- **`--tsne`**, **`--tsne-perplexity`** — tSNE of the loci over the selected samples (default perplexity 30)

## What the script does

### 1) Chromosome naming is harmonised automatically

A panel file named Ensembl-style (`1`) against a `chr1` inference BED produces **zero**
overlaps, and the failure is silent — the run completes and the column reads as an RBP that
binds nothing. Five columns were lost that way before this was fixed.

The script takes a **majority vote** over each panel file's first 2,000 data chromosomes and
rewrites only the mismatched files into its temp directory, logging:

```text
Chromosome naming: rewrote 5 panel column(s) to match the inference BED (chr-prefixed style): ...
```

Majority vote, not the first row: BEDs are ASCII-sorted, which puts `GL`/`KI` scaffolds
before `chr1`, so sampling one line misclassified 278 of 302 columns. Panel files with no
data rows are warned about, since they look identical in the output.

### 2) Build per-locus signal vectors

Each locus is expanded by `--window` (`bedtools slop`), then intersected with each panel file
**strand-aware** (`bedtools intersect -s`). For each locus the script builds a vector of
length `2*window+1` holding the summed panel **score** at each relative offset.

Offsets are strand-aligned:

- locus `+`: offset = `anchor - locus_start`
- locus `-`: offset = `-(anchor - locus_start)`, so positive offsets are always 5'→3'

### Panel anchor: `start` vs `midpoint`

A panel interval's **entire score lands on one offset**; it is not spread across its width.
`--panel-anchor` picks which offset.

For 1 nt crosslink sites the two settings are identical (`start == midpoint` when
`end == start + 1`). They diverge as soon as the panel holds **peaks**, and not harmlessly:

- a peak of width `w` centred on a `+` strand locus scores at **`-w/2`**
- the strand flip sends the same peak on a `-` strand locus to **`+w/2`**, because a peak's
  genomic start is its 3' end there

Real central binding is therefore split into a spurious **`+/-w/2` doublet with a hole at 0**.
Measured on a panel of Clippy peaks (mean width ~11 bp) against 29,018 centred loci:

| `--panel-anchor` | `+` strand median | `-` strand median | separation |
|---|---:|---:|---:|
| `start` | -5.0 | +5.0 | **10.0 nt** |
| `midpoint` | 0.0 | 0.0 | **0.0 nt** |

**Use `midpoint` for any peak-based panel.** It is a no-op for crosslink input, so it is safe
to leave on. `start` remains the default only so existing runs reproduce byte for byte.

Note `max_binding_offset` in the summary table centres a 5 nt sliding window, so it carries
its own quantisation of a couple of nt independently of this setting.

### Heatmap colour scaling

`percentile` (default) applies `log1p`, then scales against the given percentile of the
**non-zero** values and clips. Empty cells stay at 0, so the full palette carries signal.

`logistic` centres on the matrix **median**. The locus x sample matrix is sparse, so that
median is ~0 and **every empty cell maps to exactly 0.5** — mid-palette. The whole range
below 0.5 goes unused and the colourbar starts at 0.5, which is why a sparse run looked
uniformly flat no matter how the row filter was set. On a sparse test matrix (65% empty):

| scaling | range | empty cells render at | median non-zero cell |
|---|---|---:|---:|
| `logistic` | 0.50 – 1.00 | 0.50 | 0.67 |
| `percentile` (log1p, 99th) | 0.00 – 1.00 | 0.00 | 0.53 |

The `log1p` step is what makes this usable rather than merely correct. Support counts are
heavy-tailed — non-zero median 10 against a maximum of 333 on that matrix — so scaling raw
values against the 99th percentile put the median cell at **0.11**, *darker* than the
logistic scaling it replaced. After `log1p` the same settings put it at 0.53.

This also drives the **column dendrogram**, built from cosine distances on the scaled matrix:
under `logistic` every column carries a large constant 0.5 component from its empty cells,
which compresses the distances; under `percentile` they reflect actual co-occurrence. It does
**not** affect the row clusters, which come from `(matrix > 0)` binarisation.

### The heatmap row filter

Loci enter the heatmap only if their support clears `--support-pct`. Support here is the
**summed score of every overlapping panel interval, across every column and the whole
window** — so it has no intrinsic scale. It rises with panel width, with sequencing depth,
and with how bound the region is, and none of those are properties of the locus.

That is why this is a percentile rather than an absolute cut. Measured on the THRAP3 exonic
and intronic subsets, median row support was **1,594** and **296** — a 5.4x gap — so the flat
threshold of 200 used previously dropped **11.3%** of exonic loci and **43.5%** of intronic
ones. The intronic heatmap was a top-56% subset of its loci while the exonic one was a
top-89% subset, and the two were being read side by side.

Loci with **zero** support are always dropped, whatever the percentile resolves to. That
clause is load-bearing: 18.6% of the THRAP3 intronic loci have no support from any panel
column, so the 10th percentile of that set is literally 0 and a bare percentile would keep
~750 blank heatmap rows. At `--support-pct 10` the exonic set resolves to a threshold of 166
and loses 10.0%; the intronic set resolves to 0, the zero clause fires, and it loses 18.6%.
The residual asymmetry is then real emptiness rather than an artefact of scale.

Nothing is lost either way: `binf_summary.tsv` and `protein_ranking.tsv` cover every locus
and every sample, unfiltered.

### Protein selection: `total` vs `enrichment`

`--top-proteins K` picks how many panel samples reach the heatmap, the metaprofile and the
tSNE. `--protein-select` decides *how* they are picked.

`total` ranks by summed support. That is **depth-biased**: a deeply sequenced sample outranks
a shallow one whose binding is far better positioned. In a controlled test — decoys built
from a real replicate at matched depth, one shifted 40 nt and one jittered — the shifted
decoy ranked **3rd**, above genuine replicates.

`enrichment` (default) asks, per locus, whether the sample's strongest binding lands within
`--enrichment-window` of the locus, then takes the fraction of that sample's signal-bearing
loci where it does. It is depth-free.

**They are only weakly related, and neither is simply right.** On the THRAP3 run, Spearman
between the two is **+0.35** across 224 samples. The consequence is concrete: HNRNPC ranks
**1st** by enrichment but **116th** by depth, while BCLAF1 — THRAP3's known complex partner —
ranks **5th** by depth and **16th** by enrichment.

The trade-off:

- **depth** finds partners that bind abundantly nearby, but ranks a deeply sequenced sample
  above a well-positioned shallow one, and cannot tell "binds this locus" from "binds
  everything".
- **enrichment** finds partners that bind *at* the locus, but a partner binding 5–10 nt away
  scores nothing, and the metric is trivially 1.0 for a sample touching a single locus.

Because of that last point there is deliberately **no eligibility gate**. Gates on minimum
loci made two runs over different locus sets incomparable and rendered "excluded"
indistinguishable from "absent" in the output. Instead `protein_ranking.tsv` reports
`loci_with_signal` for every sample, so a top-ranked sample resting on 3 loci is visible.

Run both and compare — the ranking table is identical in shape either way, so only
`selected_for_figures` differs.

## Outputs

### `metaprofile.png` (always)

- **left axis**: `counts.mean(axis=0)`, Gaussian-smoothed. The denominator is **every** locus,
  including those where the sample has no signal, so a curve is diluted by non-binding loci
  rather than describing the sites it does bind.
- **right axis**: the same curve times the locus count, i.e. total summed panel score. One
  constant rescale, so the two axes agree pixel for pixel — the only normalisation that can
  honestly share an axis. Anything per-sample (dividing each curve by its own maximum)
  reorders the curves and needs its own panel.
- curves are the top `--top-proteins` of the `--protein-select` ranking — the **same** set the
  heatmap shows, so the two figures never disagree about which samples matter
- legend entries carry each sample's grand total; curves past the tenth switch linestyle,
  since the colour cycle is 10 long

### `binf_support_heatmap.png` (always)

- rows: the selected samples, labelled `NAME [rank]` by `--protein-select` position
- columns: loci passing `--support-pct`
- values: per-locus total support, scaled per `--heatmap-scale`
- **row order** comes from the sample dendrogram (cosine distance, average linkage), which
  groups by co-occurrence — *not* by rank. The bracketed rank makes the two orderings
  comparable, and makes it obvious when a visually dominant row is one the ranking placed
  near the cut.
- **locus order**: by total support descending; with `-n`, by cluster then support

### `protein_ranking.tsv` (always)

One row per panel sample, sorted by total peak support. Columns:

| column | meaning |
|---|---|
| `sample` | samplesheet `group` |
| `mean_peak_support` | total / n_loci — the metaprofile's left axis, integrated |
| `total_peak_support` | summed panel score across every locus and offset |
| `loci_with_signal`, `frac_loci_with_signal` | how much of the locus set this sample touches at all |
| `frac_centred` | fraction of signal-bearing loci peaking within `--enrichment-window` |
| `total_rank`, `enrichment_rank` | **both** rankings, always, whichever `--protein-select` was used |
| `selected_for_figures` | reached the top `--top-proteins` |

Because both rankings are always written, an `enrichment` run and a `total` run over the same
loci are comparable row for row — only `selected_for_figures` differs.

### `binf_summary.tsv` (always)

One row per locus, and for each sample five statistics over its `-window..+window` vector:
`_total_overlaps`, `_variance`, `_pearson_median_skew`, `_kurtosis_excess`,
`_max_binding_offset` (centre of the highest-summing 5 nt sliding window; 0 for empty loci).

The `_total_overlaps` suffix is load-bearing — `compare_binding_frequency.py`,
`plot_offset_distribution.py` and `plot_cluster_metaprofile_at_loci.py` all key off it.

### `binf_summary_tsne.png` (with `--tsne`)

Features are the `*_total_overlaps` columns for the selected samples, logistic-scaled to match
the heatmap. One point per locus. With `-n`, coloured by cluster in the **same** hues the
heatmap's cluster bar uses, so `C3` here is `C3` there; loci failing the row filter are grey.
Without `-n`, a single-colour scatter.

### `binf_heatmap_clusters.tsv` and `metaprofile_cluster_C*.png` (with `-n`)

k-means on the **binarised** matrix (`support > 0`), so clusters describe *which* samples are
present, not how much. The TSV carries `binf_chr_start_end`, `chrom`, `start`, `end`,
`row_sum_support`, `passes_heatmap_filter` and `heatmap_cluster` (`NA` if the locus failed
the filter).

This file is the **only** carrier of cluster membership, and both
`scripts/plot_offset_distribution.py --clusters` and
`scripts/plot_cluster_metaprofile_at_loci.py --clusters` read it.

One metaprofile is written per cluster, over the same samples as the global one.

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
| `run_thrap3_intersect.sbatch [enrichment\|total]` | THRAP3 against the eCLIP panel, all anchors. |
| `run_thrap3_region.sbatch <exonic\|intronic> [enrichment\|total]` | THRAP3 split by transcript region; identical settings so divergence is region, not parameters. Second argument switches the ranking and the output directory, for comparing the two. |
| `run_centro_intersect.sbatch` | Centrosome apex CLIP, uncontrolled. |
| `run_centro_controlled.sbatch <specific\|ntcontrol>` | Centrosome with APEX controls applied. Run **both**; neither is interpretable alone. |

### Python version

`intersect_inference_bed.py` needs the `rbpeek` conda env (Python 3.12). The helper scripts
above are written to run under **Python 3.6** as well, since they are typically invoked by
hand on an HPC login node where `python3` is the system interpreter.

## Notes

- Runtime is dominated by `bedtools intersect`, once per panel column.
- k-means and tSNE use a fixed seed (`RANDOM_STATE = 42`). It used to be a flag; nobody varied it, and a run that reproduces is worth more than one that can be reseeded.

