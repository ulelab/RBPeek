# `intersect_inference_bed.py`

Rank a panel of CLIP/iCLIP/eCLIP samples by how much of their binding sits at the loci of an
**inference BED**, answering "which RBPs co-occupy these sites, and with what profile shape".

Every run writes:

| file | what it is |
|---|---|
| `metaprofile.pdf` | normalised mean support vs offset, for the first 10 selected samples |
| `sample_summary.tsv` | one row per panel sample: central binding (the ranking score), rank, proportional binding, raw support and shape statistics |
| `binf_support_heatmap.pdf` | loci x selected samples, normalised support |
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

BED6+, strand in column 6; extra columns are ignored. Loci may be 1 nt anchors or intervals
such as peaks: every locus is anchored at `(start+end)//2`, which is the start itself for a
1 nt locus. Coordinates may repeat.

### Normalisation regions (`--norm-bed`)

BED6 of the regions each sample's normalising cDNA is summed over. Match it to the loci:
`regions_exonic.bed` for an exonic locus set, `regions_intronic.bed` for an intronic one.
`split_inference_bed_by_region.py` writes both.

### Genome sizes (`--genome`)

Used for `bedtools slop` when expanding loci by `--window`. Make it with `samtools faidx`
then `cut -f1,2` on the `.fa.fai`.

## CLI

```bash
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b THRAP3/THRAP3_merged_min2rep_anchors_exonic.bed \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  --norm-bed THRAP3/regions_exonic.bed \
  -o results/thrap3_exonic \
  --tsne
```

### Options

**Required**

- **`-x/--xldir`** — root the samplesheet's `file` paths resolve against
- **`-b/--bed`** — inference BED
- **`-s/--samplesheet`** — TSV with `file` and `group` columns
- **`--norm-bed`** — normalisation regions

**Counting**

- **`--window`** — half-window in nt (default 100); offsets run `-window..+window`
- **`--gaussian-sigma`** — metaprofile smoothing sigma (default 2.0)
- **`--genome`** — genome sizes file

**Selection and heatmap**

- **`--support-pct`** — keep the top P% of samples by central binding (default 30)
- **`--central-window`** — half-width in nt of the window around nt 0 counted by the ranking (default 10, i.e. −10..+10)
- **`--heatmap-scale-percentile`** — non-zero percentile mapped to the top of the colour range (default 99)

**Optional extras**

- **`-n/--n-clusters`** — k-means the heatmap's loci into this many groups; also writes `binf_heatmap_clusters.tsv` and one metaprofile per cluster, and colours the tSNE
- **`--tsne`**, **`--tsne-perplexity`** — tSNE of the heatmap's loci (default perplexity 30)

## What the script does

### 1) Per-locus signal vectors

Each locus is expanded by `--window` and intersected with each panel file **strand-aware**.
For each locus the script builds a vector of length `2*window+1` holding the summed peak
**score** at each offset. For Clippy peaks the score is the crosslink cDNA inside the peak.

Every panel interval is anchored at `(start+end)//2`, so a crosslink sits at its own position
and a peak at its centre. Offsets are strand-aligned, so positive is always 5'→3' of the locus.

### 2) Normalisation

For each sample:

- **`region_cdna`** — summed score of its peaks whose midpoint lies in `--norm-bed`, strand-aware, **chrM excluded**
- **`locus_cdna`** — summed score of its **distinct** peaks within `±window` of any locus
- **`proportional_binding`** = `locus_cdna / region_cdna`

Peaks are counted once because neighbouring loci's windows overlap,
yet `total_peak_support` still reports that sum for comparison.

### 3) Selection

Samples are ranked by **`central_binding`**: for each locus, the support at offsets within
±`--central-window` of nt 0 (default ±10, 21 nt) is summed, divided by the sample's `region_cdna`
and multiplied by 10⁶, then `log1p`-transformed. The score is the mean of that over **all** loci. The top `--support-pct`%
(68 of 224 at the default 30) go to the heatmap, tSNE and clustering, and the metaprofile draws
the first 10 of them.

- **Central window** — only binding at the locus counts; binding further out does not affect the rank.
- **÷ region cDNA** — sequencing depth nromalisation
- **log1p** — no handful of loci can decide the rank, this is ranked by raw mean support.
- **mean over all loci** — loci a sample doesn't bind score 0, which keeps sparse samples down.

A drawback is a sample binding a few inference loci very strongly ranks
below one binding many sites moderately. `proportional_binding` and `mean_peak_support` are
still reported for every sample.
Figures show per-locus support ÷ `region_cdna` × 10⁶, labelled "per M region cDNA".

Loci with no support from any selected sample are dropped from the heatmap; every other locus
is kept.

### Heatmap colour scaling

The matrix is `log1p`-transformed, then scaled against the `--heatmap-scale-percentile` of its
**non-zero** values and clipped. Empty cells stay at 0, so the whole palette carries signal.
The percentile is over non-zero values because on a sparse matrix a percentile over all cells
is itself 0. The `log1p` matters: support is heavy-tailed, and scaling raw values against the
99th percentile put the median cell at 0.11.

## Outputs

### `metaprofile.pdf`

- **left axis**: mean normalised support per locus. Every locus is in the denominator, including
  those where the sample has no signal.
- **right axis**: the same curve times the locus count. One constant, so both axes describe the
  same pixels.
- curves are the first 10 selected samples; legend entries give each one's central
  binding. Curves past the tenth switch linestyle, since the colour cycle is 10 long.

### `binf_support_heatmap.pdf`

- title: the inference BED name, how many samples were selected out of how many, and the
  number of loci with support
- rows: the selected samples, labelled `NAME [rank]`
- columns: loci with support from at least one selected sample
- values: normalised support, colour-scaled as above
- **row order** comes from the sample dendrogram (cosine distance, average linkage), which
  groups by co-occurrence, not rank — hence the bracketed rank on each label
- **locus order**: by total support; with `-n`, by cluster then support

### `sample_summary.tsv`

One row per panel sample, sorted by rank.

| column | meaning |
|---|---|
| `sample` | samplesheet `group` |
| `rank`, `selected` | rank by `central_binding`; whether it made the top `--support-pct`% |
| `central_binding` | the ranking score; see [Selection](#3-selection) |
| `proportional_binding`, `locus_cdna`, `region_cdna` | see [Normalisation](#2-normalisation) |
| `mean_peak_support`, `total_peak_support` | raw per-window support: total ÷ number of loci, and the sum over every locus and offset |
| `loci_with_signal`, `frac_loci_with_signal` | how much of the locus set the sample touches at all |
| `mean_binding_offset` | mean over signal-bearing loci of the centre of the strongest 5 nt window |
| `mean_variance` | mean per-locus variance, on the normalised scale so samples of different depth compare |
| `mean_pearson_skew`, `mean_kurtosis` | means over signal-bearing loci |

### `binf_summary_tsne.png` (with `--tsne`)

One point per heatmap locus, on the same scaled matrix the heatmap draws. With `-n`, coloured
by cluster in the same hues as the heatmap's cluster bar.

### `binf_heatmap_clusters.tsv` and `metaprofile_cluster_C*.pdf` (with `-n`)

k-means on the **binarised** matrix (`support > 0`), so clusters describe *which* samples are
present, not how much. The TSV carries `binf_chr_start_end`, `chrom`, `start`, `end`,
`row_sum_support`, `passes_heatmap_filter` and `heatmap_cluster` (`NA` for loci not in the
heatmap). One metaprofile is written per cluster, over the same samples as the global one.

## Repository layout

`intersect_inference_bed.py` is the analysis engine. Everything else prepares its inputs.

| script | purpose |
|---|---|
| `build_inference_bed_from_peaks.py` | **General.** Replicate Clippy peak calls -> inference BED. Normalises chromosome names, merges strand-aware, keeps regions supported by `--min-reps` replicates, collapses each to a 1 nt midpoint anchor. `--subtract` drops regions overlapping control peaks, for proximity-labelling baits. |
| `build_thrap3_inference_bed.py` | The THRAP3-specific original, kept for reproducibility of that analysis. Superseded by the general script above. |
| `merge_replicate_peaks.py` | Pools replicate peak BEDs by SUMMING scores at identical chrom/start/end/strand, then reformats to BED6. Exact-match grouping, so it pools signal rather than assessing reproducibility - for the latter use the overlap-based merge above. |
| `split_inference_bed_by_region.py` | Splits an inference BED into **exonic** and **intronic** subsets from a GTF (strand-aware, exon-priority, so the two are disjoint), and writes `regions_exonic.bed` / `regions_intronic.bed` for `--norm-bed`. `--drop-chrM` removes mitochondrial anchors first. Defaults to GENCODE v39. |
| `run_thrap3_region.sbatch <exonic\|intronic>` | THRAP3 against the eCLIP panel for one region. |

### Python version

`intersect_inference_bed.py` needs the `rbpeek` conda env (Python 3.12). The other scripts
also run under **Python 3.6**, since they are typically invoked by hand on an HPC login node
where `python3` is the system interpreter.

## Notes

- Figures are PDFs with embedded TrueType fonts, so text stays editable. The heatmap's cell mesh
  is rasterised at 200 dpi inside the PDF; everything else is vector.
- Runtime is dominated by `bedtools intersect`, twice per panel column (locus windows and normalisation regions).
- k-means and tSNE use a fixed seed (`RANDOM_STATE = 42`), so runs reproduce.
