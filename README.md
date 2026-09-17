# RBPeek

Ranks a panel of CLIP samples by how much of their binding sits at the loci of an **inference
BED** (here, reproducible THRAP3 sites), then plots the top-ranked samples.

## Quick start (THRAP3)

Once, split the anchors into exonic and intronic sets and write the normalisation regions:

```bash
python3 scripts/split_inference_bed_by_region.py -b THRAP3/THRAP3_merged_min2rep_anchors.bed -o THRAP3 --drop-chrM
```

Then run each region; results go to `results/thrap3_<region>/`:

```bash
sbatch --job-name=thrap3_exonic scripts/run_thrap3_region.sbatch exonic
```

```bash
sbatch --job-name=thrap3_intronic scripts/run_thrap3_region.sbatch intronic
```

The engine needs the `rbpeek` conda env (`rbpeek.yml`, Python 3.12). The other scripts also run
under Python 3.6 on the login node.

## Outputs

| file | contents |
|---|---|
| `sample_summary.tsv` | one row per panel sample, sorted by rank ([columns](#sample_summarytsv)) |
| `metaprofile.pdf` | mean normalised profile of each locus's strongest central peak vs offset, first 10 selected samples |
| `binf_support_heatmap.pdf` | loci × selected samples, each cell the strongest central peak; titled with the BED name, selection and locus count |
| `binf_summary_tsne.png` | tSNE of the heatmap loci (`--tsne`) |
| `binf_heatmap_clusters.tsv`, `metaprofile_cluster_C*.pdf` | k-means clusters of the heatmap loci, one metaprofile per cluster (`-n`) |

## Inputs and options

| flag | default | meaning |
|---|---|---|
| `-x/--xldir` | required | root the samplesheet's `file` paths resolve against |
| `-s/--samplesheet` | required | TSV of `file` and `group` (label); THRAP3 runs use `THRAP3/RBPeekSamplesheet_eCLIP.tsv` |
| `-b/--bed` | required | inference BED6+, strand in column 6; each locus is anchored at `(start+end)//2` |
| `--norm-bed` | required | regions for each sample's normalising cDNA; use the `regions_<region>.bed` matching the loci |
| `--genome` | HPC hg38 | chromosome sizes for `bedtools slop` |
| `-o/--outdir` | `results` | output directory |
| `--window` | 100 | half-window (nt) around each locus |
| `--central-window` | 10 | half-width (nt) counted by the ranking score |
| `--support-pct` | 30 | top P% of samples sent to the heatmap, tSNE and clustering |
| `--gaussian-sigma` | 2 | metaprofile smoothing |
| `--heatmap-scale-percentile` | 99 | non-zero percentile mapped to the top of the colour scale |
| `-n/--n-clusters` | off | k-means the heatmap loci on presence/absence |
| `--tsne`, `--tsne-perplexity` | off, 30 | tSNE of the heatmap loci |

Chromosome naming (`chr1` vs `1`) is harmonised automatically between the inference BED and the
panel.

## Method

1. **Signal vectors.** Each locus is expanded to ±`--window` and intersected strand-aware with
   each panel file. A peak's cDNA (Clippy score) is spread evenly across its width, so a peak
   partly inside a window counts only for the overlap: an 11 nt peak centred at +12 adds 4/11 of
   its cDNA to offsets +7..+10. A 1 nt crosslink keeps its whole score at its own position.
   Offsets are strand-aligned; positive is downstream of the locus.
2. **Normalisation.**
   - `region_cdna`: the sample's peak cDNA inside `--norm-bed`, strand-aware, chrM excluded.
     Each peak counts by the fraction of its width inside the regions, which are merged first.
   - `locus_cdna`: its cDNA inside any locus window, each distinct peak counted once.
   - `proportional_binding` = `locus_cdna / region_cdna`.
3. **Ranking.** `central_binding` = mean over **all** loci of
   `log1p(support within ±central-window × 10⁶ / region_cdna)`. Every peak in the window counts.
   - The central window means only binding at the locus counts.
   - Dividing by `region_cdna` removes sequencing depth.
   - `log1p` stops a few very strong loci from deciding the rank.
   - Averaging over all loci scores unbound loci as 0. That keeps sparse samples down, but a
     sample binding a few loci very strongly ranks below one binding many loci moderately.

   Samples with no region cDNA or no support at the loci are never selected.
4. **Figures.** These draw one match per locus; the ranking above still counts every peak.
   - At each locus, only the sample's **strongest central peak** is drawn: the peak with the most
     cDNA inside ±central-window. Values are × 10⁶ / `region_cdna`.
   - A heatmap cell is that peak's cDNA inside the window. The metaprofile averages that peak
     alone, spread over its width, so it can extend past ±central-window but other peaks at the
     locus are not drawn.
   - The heatmap is `log1p`-transformed, scaled to the `--heatmap-scale-percentile` of non-zero
     cells, and drops loci with no central peak from any selected sample.
   - Heatmap rows follow a cosine-distance sample dendrogram, and the bracketed number is the
     rank. Loci are ordered by summed cell value, or by cluster with `-n`.
   - The metaprofile's right axis is its left axis × the number of loci.
   - PDFs embed TrueType fonts, and heatmap cells are rasterised at 200 dpi. k-means and tSNE use
     a fixed seed.

### `sample_summary.tsv`

| column | meaning |
|---|---|
| `rank`, `selected` | rank by `central_binding`; whether it was selected |
| `central_binding` | the ranking score |
| `proportional_binding`, `locus_cdna`, `region_cdna` | see Method 2 |
| `total_peak_support`, `mean_peak_support` | support summed over every locus window (a peak near two loci counts twice), and ÷ number of loci |
| `loci_with_signal`, `frac_loci_with_signal` | loci with any support |
| `mean_binding_offset`, `mean_variance`, `mean_pearson_skew`, `mean_kurtosis` | means over loci with signal; offset is the centre of the strongest 5 nt window; variance is on the normalised scale |

## THRAP3 inference loci

The inputs are four HEK293 HA-THRAP3 iCLIP replicates (Flow project `788995297969977723`,
CLIP-Seq v1.7, GRCh38) in `THRAP3/raw/`. `build_thrap3_inference_bed.py`:
- converts chromosome names
- merges overlapping peaks strand-aware across replicates
- keeps regions supported by ≥2 of 4 replicates
- collapses each region to its 1 nt midpoint

That gives 29,018 loci in `THRAP3_merged_min2rep_anchors.bed`, with score = replicate support.

**Caveats.** The panel is HepG2/K562 eCLIP, so cell line is confounded with RBP: compare ranks,
not absolute values. The eCLIP sheet excludes `HNRNPC-Hela-iCLIP`, the only assay-matched sample.
`RBPeekSamplesheet.tsv` has all 299 samples.

## Scripts

| script | purpose |
|---|---|
| `intersect_inference_bed.py` | the engine |
| `split_inference_bed_by_region.py` | exonic/intronic split (strand-aware, exon-priority) and the normalisation region BEDs; `--drop-chrM` |
| `run_thrap3_region.sbatch` | runs the engine for one region |
| `build_thrap3_inference_bed.py` | THRAP3 replicate peaks → inference BED |
| `merge_replicate_peaks.py` | pools replicate peak BEDs by summing scores at identical intervals (the `peaks/merged/` panel entries) |

`Centrosome/` and `decoys/` hold data from earlier analyses; the current method doesn't use them.
