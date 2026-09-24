# RBPeek

RBPeek quantifies the co-occupancy of RNA-binding proteins (RBPs) at a defined set of genomic
loci. Given an **inference BED** (here, reproducible THRAP3 binding sites) and a panel of CLIP
peak files, it ranks the panel samples by their depth-normalised binding at those loci and
plots the binding profiles of the highest-ranked samples.

## Requirements

`intersect_inference_bed.py` requires Python ≥ 3.9, `bedtools`, `numpy`, `scipy`, `matplotlib`,
`seaborn` and `scikit-learn`; `rbpeek.yml` defines a conda environment that provides them. The
other scripts require only Python ≥ 3.6 and `bedtools`.

## Usage

1. Partition the inference loci into exonic and intronic sets and generate the normalisation
   regions:

   ```bash
   python3 scripts/split_inference_bed_by_region.py -b THRAP3/THRAP3_merged_min2rep_anchors.bed -g <annotation.gtf> -o THRAP3 --drop-chrM
   ```

2. Run the analysis for one region:

   ```bash
   python3 scripts/intersect_inference_bed.py \
     -x <panel directory> \
     -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
     -b THRAP3/THRAP3_merged_min2rep_anchors_exonic.bed \
     --norm-bed THRAP3/regions_exonic.bed \
     --genome <chromosome sizes> \
     -o results/thrap3_exonic \
     --support-pct 20 --tsne
   ```

   `scripts/run_thrap3_region.sbatch <exonic|intronic>` runs this step under SLURM with the
   settings used for the THRAP3 analysis.

## Options

| option | default | description |
|---|---|---|
| `-x/--xldir` | required | directory against which the samplesheet's `file` paths are resolved |
| `-s/--samplesheet` | required | TSV with columns `file` and `group` (sample label) |
| `-b/--bed` | required | inference BED6+, strand in column 6; each locus is anchored at `(start + end) // 2` |
| `--norm-bed` | primary chromosomes of `--genome`, both strands | BED6 of the regions over which each sample's normalising cDNA is summed (`regions_exonic.bed` for exonic loci, `regions_intronic.bed` for intronic loci) |
| `--genome` | required | chromosome sizes file for `bedtools slop` |
| `-o/--outdir` | `results` | output directory |
| `--window` | 100 | half-width (nt) of the window around each locus |
| `--central-window` | 10 | half-width (nt) of the central window used for ranking |
| `--support-pct` | 30 | percentage of top-ranked samples shown in the heatmap and tSNE |
| `--highlight` | none | sample labels, or substrings of them, drawn in black on the metaprofile in addition to the ten highest-ranked samples |
| `--heatmap-scale-percentile` | 99 | percentile of non-zero heatmap values mapped to the top of the colour scale |
| `--tsne`, `--tsne-perplexity` | off, 30 | tSNE embedding of the heatmap loci |

Chromosome naming (`chr1` or `1`) is harmonised automatically between the inference BED and
the panel files.

## Method

1. **Signal.** Each locus is extended by `--window` nt on either side and intersected with each
   panel file on the same strand. The cDNA count of a peak (BED score) is distributed uniformly
   across its width, so a peak that partly overlaps a window contributes only the overlapping
   fraction, and a 1 nt crosslink site retains its full score at its own position. Offsets are
   strand-aligned, with positive values downstream of the locus. This yields, for each sample,
   a locus × offset support matrix *c*(*l*, *o*).
2. **Normalisation.** `region_cdna` is the cDNA of the sample's peaks within `--norm-bed` on
   the same strand (by default, the primary chromosomes of `--genome` on both strands). Each peak is weighted by the fraction of its width inside the regions,
   overlapping regions are merged beforehand, and mitochondrial peaks are excluded because
   mitochondrial rRNA is a major source of background in eCLIP libraries. Normalised support
   is *x*(*l*, *o*) = *c*(*l*, *o*) × 10⁶ / `region_cdna`.
3. **Metaprofile.** For each sample, *m*(*o*) is the mean over all loci of
   log(1 + *x*(*l*, *o*)). The logarithm is applied per locus before averaging, so that the
   curve reflects the breadth of binding across loci and is not dominated by a small number of
   strongly bound loci. Loci without signal contribute 0.
4. **Ranking.** `central_binding` is the area under *m*(*o*) within ±`--central-window` nt.
   Samples are ranked by `central_binding`, and the top `--support-pct` percent are selected.
   Because the ranking score and the plotted curve are the same statistic, curve height reflects
   rank. Samples without region cDNA, or without a peak in the central window of any locus, are
   not selected.
5. **Figures.**
   - `metaprofile.pdf` shows *m*(*o*) for the ten highest-ranked samples, without smoothing, so
     the area between the red dotted lines equals `central_binding`. These lines mark the central
     window, and the legend gives the rank and `central_binding` of each sample. Samples named
     with `--highlight` are added to the plot and drawn in black.
   - `binf_support_heatmap.pdf` shows loci × selected samples. Each cell is the strongest single
     peak of the sample within the central window of the locus (the cDNA of that peak inside the
     window, per million region cDNA), log(1 + *x*)-transformed and scaled to the
     `--heatmap-scale-percentile` of the non-zero cells. Samples are ordered by hierarchical
     clustering (cosine distance, average linkage) and labelled with their rank; loci are
     ordered by summed cell value, and loci without a central peak in any selected sample are
     omitted.
   - `binf_summary_tsne.png` is a tSNE embedding of the heatmap loci (fixed random seed).

The width of the central peak in the metaprofile is determined by the panel peak calls and
should not be interpreted as a binding footprint: the panel peaks are approximately 11 nt wide,
and the cDNA of each peak is distributed uniformly across its width.

## Outputs

| file | contents |
|---|---|
| `metaprofile.pdf` | *m*(*o*) of the ten highest-ranked selected samples |
| `sample_summary.tsv` | one row per panel sample, sorted by rank |
| `binf_support_heatmap.pdf` | loci × selected samples |
| `binf_summary_tsne.png` | tSNE embedding of the heatmap loci (`--tsne`) |

Columns of `sample_summary.tsv`:

| column | description |
|---|---|
| `sample` | sample label from the samplesheet |
| `rank`, `selected` | rank by `central_binding`; whether the sample was selected |
| `central_binding` | area under *m*(*o*) within the central window |
| `region_cdna` | normalising cDNA (Method 2) |
| `locus_cdna`, `proportional_binding` | cDNA of the distinct peaks within any locus window, each weighted by the fraction of its width inside; and its ratio to `region_cdna` |
| `total_peak_support`, `mean_peak_support` | support summed over all locus windows (a peak within the windows of two loci is counted in both), and that sum divided by the number of loci |
| `loci_with_signal`, `frac_loci_with_signal` | number and fraction of loci with any support within the window |
| `mean_binding_offset` | mean, over loci with signal, of the centre of the 5 nt window with the highest support |
| `mean_variance`, `mean_pearson_skew`, `mean_kurtosis` | means, over loci with signal, of the variance (normalised scale), Pearson median skewness and excess kurtosis of the support vector |

## THRAP3 data

The inference loci derive from four HEK293 HA-THRAP3 iCLIP replicates (GRCh38). Crosslink sites
were derived from read 1 of each pair, and the Clippy peak calls are stored in `THRAP3/raw/`. `build_thrap3_inference_bed.py` converts chromosome
names to UCSC style, merges overlapping peaks across replicates on each strand, retains regions
supported by at least two of the four replicates, and reduces each region to its 1 nt midpoint.
This yields 26,030 loci (`THRAP3_merged_min2rep_anchors.bed`; the score column records
replicate support), of which 216 are mitochondrial and are removed before the loci are
partitioned into exonic and intronic sets.

The panel (`THRAP3/RBPeekSamplesheet_eCLIP.tsv`) comprises 224 HepG2 and K562 eCLIP samples;
`RBPeekSamplesheet.tsv` additionally lists the iCLIP and PAR-CLIP samples. Panel entries under
`peaks/merged/` were produced with `merge_replicate_peaks.py`.

**Limitations.** THRAP3 was profiled in HEK293 cells, whereas the panel consists of HepG2 and
K562 data, so cell line is confounded with RBP identity; ranks should be compared, not
absolute values. `central_binding` increases with the fraction of loci at which a sample has
signal, and the number of peaks called depends on library depth, so normalisation by
`region_cdna` does not remove the influence of depth entirely.

## Centrosome data

`build_centrosome_inference_bed.py` builds inference loci for proximity-CLIP (APEX) data, in
which peak cDNA reflects the proportion of a transcript in the labelled compartment rather
than a discrete binding site. Peak cDNA is assigned to gene bodies (each peak weighted by the
fraction of its width inside the gene, same strand) and expressed per million of the sample's
non-mitochondrial peak cDNA. Genes are retained if they pass the differential-expression
filter (`log2FoldChange > 1`, `pvalue <= 0.05`) and the ratio of the experimental to the
control group mean, with a pseudocount of 1 per million added to both, is at least 2. The
anchors are the midpoints of the pooled experimental peaks that lie within retained genes.

```bash
python3 scripts/build_centrosome_inference_bed.py --deseq <DESeq2 table.xlsx> -g <annotation.gtf>
```

Because the loci are not partitioned by region, `intersect_inference_bed.py` is run without
`--norm-bed`, so that each sample is normalised by its peak cDNA on the primary chromosomes.

## Scripts

| script | purpose |
|---|---|
| `intersect_inference_bed.py` | the analysis described above |
| `split_inference_bed_by_region.py` | exonic and intronic locus sets (strand-aware, exons given priority) and the corresponding normalisation regions; `--drop-chrM` removes mitochondrial loci |
| `build_thrap3_inference_bed.py` | replicate THRAP3 peak calls → inference BED |
| `build_centrosome_inference_bed.py` | proximity-CLIP peak calls and differential-expression table → inference BED of enriched genes |
| `merge_replicate_peaks.py` | pools replicate peak BEDs by summing scores at identical intervals |
| `run_thrap3_region.sbatch` | SLURM job for one region of the THRAP3 analysis |
