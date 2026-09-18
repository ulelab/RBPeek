# RBPeek

RBPeek quantifies the co-occupancy of RNA-binding proteins at a defined set of genomic loci. Given
an **inference BED** (here, reproducible THRAP3 binding sites) and a panel of CLIP peak files, it
ranks the panel samples by their depth-normalised binding at those loci and plots the binding
profiles of the highest-ranked samples.

## Quick start (THRAP3)

Partition the anchors into exonic and intronic sets and generate the normalisation regions (run
once):

```bash
python3 scripts/split_inference_bed_by_region.py -b THRAP3/THRAP3_merged_min2rep_anchors.bed -o THRAP3 --drop-chrM
```

Run each region; results are written to `results/thrap3_<region>/`:

```bash
sbatch --job-name=thrap3_exonic scripts/run_thrap3_region.sbatch exonic
```

```bash
sbatch --job-name=thrap3_intronic scripts/run_thrap3_region.sbatch intronic
```

`intersect_inference_bed.py` requires the `rbpeek` conda environment (`rbpeek.yml`, Python 3.12)
and `bedtools`. The split and build scripts are also compatible with Python 3.6.

## Outputs

| file | contents |
|---|---|
| `sample_summary.tsv` | one row per panel sample, sorted by rank ([columns](#sample_summarytsv)) |
| `metaprofile.pdf` | mean log-transformed, depth-normalised support as a function of offset across ±`--window`, for the 10 highest-ranked selected samples |
| `metaprofile_maxpeak.pdf`, `metaprofile_windowfrac.pdf` | shape-only views of the same samples (Method 4) |
| `metaprofile_matrix.npz` | raw per-locus, per-offset support of the top-ranked samples (`--save-matrix`) |
| `binf_support_heatmap.pdf` | loci × selected samples; each cell is the strongest central peak. The title gives the BED name, the selection and the locus count |
| `binf_summary_tsne.png` | tSNE embedding of the heatmap loci (`--tsne`) |
| `binf_heatmap_clusters.tsv`, `metaprofile_cluster_C*.pdf` | k-means clusters of the heatmap loci and one metaprofile per cluster (`-n`) |

## Inputs and options

| flag | default | meaning |
|---|---|---|
| `-x/--xldir` | required | root directory against which the samplesheet's `file` paths are resolved |
| `-s/--samplesheet` | required | TSV with columns `file` and `group` (sample label); the THRAP3 runs use `THRAP3/RBPeekSamplesheet_eCLIP.tsv` |
| `-b/--bed` | required | inference BED6+, strand in column 6; each locus is anchored at `(start+end)//2` |
| `--norm-bed` | required | regions over which each sample's normalising cDNA is summed; use the `regions_<region>.bed` matching the loci |
| `--genome` | HPC hg38 | chromosome sizes for `bedtools slop` |
| `-o/--outdir` | `results` | output directory |
| `--window` | 100 | half-window (nt) around each locus |
| `--central-window` | 10 | half-width (nt) of the window scored for ranking |
| `--support-pct` | 30 | percentage of top-ranked samples passed to the heatmap, tSNE and clustering (the THRAP3 runner uses 20) |
| `--gaussian-sigma` | 2 | standard deviation (nt) of the metaprofile smoothing kernel |
| `--heatmap-scale-percentile` | 99 | percentile of non-zero cells mapped to the top of the colour scale |
| `--save-matrix [N]` | off (N = 10) | write the raw counts of the N top-ranked samples to `metaprofile_matrix.npz` |
| `-n/--n-clusters` | off | number of k-means clusters of the heatmap loci (presence/absence) |
| `--tsne`, `--tsne-perplexity` | off, 30 | tSNE of the heatmap loci |

Chromosome naming (`chr1` vs `1`) is harmonised automatically between the inference BED and the
panel files.

## Method

1. **Signal vectors.** Each locus is extended to ±`--window` and intersected, strand-aware, with
   each panel file. The cDNA count of a peak (Clippy score) is distributed uniformly across its
   width, so a peak that partly overlaps a window contributes only the overlapping fraction: an
   11 nt peak centred at +12 contributes 4/11 of its cDNA to offsets +7 to +10. A 1 nt crosslink
   site retains its full score at its own position. Offsets are strand-aligned, with positive
   values downstream of the locus.
2. **Normalisation.**
   - `region_cdna`: the sample's peak cDNA within `--norm-bed` (strand-aware, chrM excluded).
     Each peak is weighted by the fraction of its width inside the regions, which are merged
     beforehand.
   - `locus_cdna`: the sample's cDNA within any locus window, with each distinct peak counted
     once.
   - `proportional_binding` = `locus_cdna / region_cdna`.
3. **Ranking.** `central_binding` is the mean over **all** loci of
   `log1p(support within ±central-window × 10⁶ / region_cdna)`, where every peak in the window
   contributes.
   - Restricting the score to the central window limits it to binding at the locus itself.
   - Division by `region_cdna` normalises for sequencing depth.
   - The `log1p` transform limits the influence of a small number of very strongly bound loci.
   - Averaging over all loci assigns unbound loci a value of 0. This penalises sparse samples;
     consequently, a sample that binds few loci strongly ranks below one that binds many loci
     moderately.

   Samples with no region cDNA, or with no peak inside the central window of any locus, are
   excluded from selection.
4. **Figures.**
   - The metaprofile shows, at each offset across ±`--window`, the mean over all loci of
     `log1p(support × 10⁶ / region_cdna)`, every peak included, smoothed with a Gaussian kernel.
     This applies the ranking's transformation per offset, so curve height follows rank. A
     linear mean of normalised support does not, because a few very strongly bound loci set
     its height, whereas the ranking rewards binding many loci. Red dotted lines mark
     ±central-window, the legend gives each sample's rank and score, and the title names the
     inference BED (and therefore the region). The plot area is square with 12 pt text.
   - Two shape-only views are also written. Both give every curve the same overall size, so
     height carries no information about rank; they compare how sharply binding is centred.
     - **maxpeak**: the metaprofile curve divided by its own maximum.
     - **windowfrac**: each locus's profile is divided by its own total over ±`--window`, then
       averaged over the loci the sample binds. Every bound locus counts equally, and depth and
       `region_cdna` cancel. The dashed line marks no positional preference (1 / window width).
   - `--save-matrix` stores the raw counts (no normalisation, no log) as sparse triplets, so
     normalisations can be tried without rerunning:
     `meta, counts = load_counts_matrix("metaprofile_matrix.npz")` (in
     `intersect_inference_bed.py`) returns the sample names, ranks, `central_binding`,
     `region_cdna`, offsets and locus names, and one loci × offsets array per sample.
   - Heatmap values are scaled by 10⁶ / `region_cdna`.
   - Each heatmap cell is the sample's **strongest central peak** at that locus: the cDNA inside
     ±central-window of the peak contributing the most cDNA there. The ranking, in contrast,
     sums all peaks.
   - The heatmap is `log1p`-transformed and scaled to the `--heatmap-scale-percentile` of
     non-zero cells. Loci without a central peak from any selected sample are omitted.
   - Heatmap rows are ordered by hierarchical clustering of samples (cosine distance, average
     linkage); the bracketed number is the sample's rank. Loci are ordered by summed cell value,
     or by cluster when `-n` is given.
   - PDFs embed TrueType fonts, and the heatmap cells are rasterised at 200 dpi. k-means and
     tSNE use a fixed random seed.

### `sample_summary.tsv`

| column | meaning |
|---|---|
| `rank`, `selected` | rank by `central_binding`; whether the sample was selected |
| `central_binding` | the ranking score (Method 3) |
| `proportional_binding`, `locus_cdna`, `region_cdna` | see Method 2 |
| `total_peak_support`, `mean_peak_support` | support summed over every locus window (a peak near two loci is counted in both), and that sum divided by the number of loci |
| `loci_with_signal`, `frac_loci_with_signal` | loci with any support within ±`--window` |
| `mean_binding_offset`, `mean_variance`, `mean_pearson_skew`, `mean_kurtosis` | means over loci with signal; the offset is the centre of the strongest 5 nt window, and the variance is on the normalised scale |

## THRAP3 inference loci

The inference loci derive from four HEK293 HA-THRAP3 iCLIP replicates (Flow project
`788995297969977723`, CLIP-Seq v1.7, GRCh38), stored in `THRAP3/raw/`.
`build_thrap3_inference_bed.py`:
- converts chromosome names to UCSC style
- merges overlapping peaks across replicates, strand-aware
- retains regions supported by at least 2 of the 4 replicates
- reduces each region to its 1 nt midpoint

This yields 29,018 loci in `THRAP3_merged_min2rep_anchors.bed`, with the score column recording
replicate support. With `--drop-chrM`, the split assigns 19,917 loci to exons and 8,666 to
introns.

**Limitations.** The panel consists of HepG2 and K562 eCLIP data, whereas THRAP3 was profiled in
HEK293, so cell line is confounded with RBP identity; ranks should be compared, not absolute
values. The eCLIP samplesheet excludes `HNRNPC-Hela-iCLIP`, the only assay-matched sample;
`RBPeekSamplesheet.tsv` lists all 299 samples.

## Scripts

| script | purpose |
|---|---|
| `intersect_inference_bed.py` | the analysis described above |
| `split_inference_bed_by_region.py` | exonic/intronic partition of the loci (strand-aware, exon priority) and the normalisation region BEDs; `--drop-chrM` removes mitochondrial loci |
| `run_thrap3_region.sbatch` | runs the analysis for one region |
| `build_thrap3_inference_bed.py` | THRAP3 replicate peaks → inference BED |
| `merge_replicate_peaks.py` | pools replicate peak BEDs by summing scores at identical intervals (used for the `peaks/merged/` panel entries) |
| `build_gene_matrix_from_summaries.py` | Flow `*.summary_gene.tsv` files → gene × sample count and RPKM matrices |
| `plot_expression_heatmap.py` | clustered expression heatmap for a gene list (`.xlsx` or text); plots the 100 most variable matched genes by default |

`Centrosome/` and `decoys/` contain data from earlier analyses and are not used by the THRAP3
workflow. The Centrosome gene-expression heatmap is generated with the last two scripts from the
files in `Centrosome/gene_counts/`.
