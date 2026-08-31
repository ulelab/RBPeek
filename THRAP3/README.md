# THRAP3 vs the RBP panel

Comparing reproducible THRAP3 binding sites against the 297-sample Clippy peak panel
under `../CLIP`, to ask which RBPs co-occupy THRAP3 sites and with what profile shape.

## Source data

Flow project [HA_THRAP3_CLIP](https://app.flow.bio/projects/788995297969977723)
(`788995297969977723`), "HA pulldown of induced THRAP3" — Karen Davey / Jernej Ule, DRI.

- Execution `188723932105002059`, pipeline **CLIP-Seq v1.7**, fileset **GRCh38**, status OK.
- Four iCLIP samples, all `condition=induced`, all **HEK293**, purified with `anti HA cs3724`.
  `THRAP3_L` differs from the others only in barcode and RT primer, so it is a genuine
  fourth replicate rather than a separate condition.

| Flow sample | Replicate | Peaks | Genome peak file |
|---|---|---:|---|
| THRAP3_1 | R1 | 25,139 | `THRAP3_R1_genome.clippy._..._Peaks.bed` |
| THRAP3_2 | R2 | 9,455 | `THRAP3_R2_genome.clippy._..._Peaks.bed` |
| THRAP3_3 | R3 | 49,477 | `THRAP3_R3_genome.clippy._..._Peaks.bed` |
| THRAP3_L | R4 | 45,991 | `THRAP3_R4_genome.clippy._..._Peaks.bed` |

Downloaded via `GET /api/downloads/<data_id>/<filename>` into `raw/`. The matching
`_Summits.bed` files were *not* used: at 350–550 MB each they are full rolling-mean
summit tracks, not peak summits.

## Building the inference BED

```bash
python3 scripts/build_thrap3_inference_bed.py
```

Three transformations are required before these peaks can serve as an inference BED:

1. **Chromosome naming.** Flow ran Clippy against `Homo_sapiens.GRCh38.fasta.fai`, so peaks
   are Ensembl-style (`1`). The panel and `decoys/` are UCSC-style (`chr1`).
   `intersect_inference_bed.py` normalises names only on its internal merge path, never for
   the inference BED — so without this step the run would silently intersect nothing.
2. **Reproducibility.** The libraries differ ~5x in depth, so a union is dominated by
   whichever was sequenced deepest. Peaks are merged strand-aware and kept at **>=2/4**
   replicates.
3. **1 nt anchors.** `intersect_inference_bed.py:356` computes
   `offset = xl_start - binf_site_start`, so loci wider than 1 bp smear the metaprofile.
   Each merged region collapses to its midpoint.

### Why >=2/4

Pairwise recovery (% of row's peaks overlapped by column, strand-aware):

|  | ->R1 | ->R2 | ->R3 | ->R4 |
|---|---|---|---|---|
| **R1** (25,139) | – | 24.7% | 61.0% | 58.1% |
| **R2** (9,455) | 65.3% | – | 77.9% | 72.3% |
| **R3** (49,477) | 31.1% | 15.1% | – | 42.0% |
| **R4** (45,991) | 31.9% | 15.0% | 45.1% | – |

R2 is shallow but high-precision — 65–78% of its peaks are recovered by the others, so it
is underpowered rather than noisy. R3/R4 carry many depth-driven private peaks. A union
would inherit R3/R4's depth bias; `>=3/4` would discard genuine signal and penalise R2's
contribution. `>=2/4` is the balance point.

Replicate support across 80,471 merged regions:

| support | regions | cumulative |
|---|---:|---:|
| 1/4 | 51,453 | 100% |
| 2/4 | 15,544 | 36.1% |
| 3/4 | 8,367 | 16.7% |
| 4/4 | 5,107 | 6.3% |

**Output:** `THRAP3_merged_min2rep_anchors.bed` — 29,018 loci, BED6, 1 nt,
`score` = replicate support (2/3/4), so results can be split by reproducibility tier.
Merged regions max 65 bp wide, so no over-chaining. Comparable in scale to the
15,172-row decoys BED.

## Running

```bash
sbatch scripts/run_thrap3_intersect.sbatch
```

THRAP3 is the **inference BED (rows)**, not a panel column — it is absent from
`RBPeekSamplesheet.tsv`, so there is no self-comparison to exclude. Output to
`results/thrap3/`.

## Priority RBPs

| RBP | Panel columns |
|---|---|
| HNRNPC | `HepG2-HNRNPC-eCLIP`, `K562-HNRNPC-eCLIP`, `HNRNPC-Hela-iCLIP`, `HNRNPC-PARCLIP` |
| ELAVL1 | `ELAVL1-PARCLIP` |
| SRSF7 | `HepG2-SRSF7-eCLIP`, `K562-SRSF7-eCLIP` |

HNRNPC appears in four independent assay/cell combinations — if its co-binding signal is
real it should be consistent across all four, which makes it a free internal control.

**The runs now use an eCLIP-only panel.** `THRAP3/RBPeekSamplesheet_eCLIP.tsv` keeps the 224
`-eCLIP` columns of 299, dropping 68 PAR-CLIP and 7 iCLIP. Assay is then uniform, so a
difference between two columns cannot be an assay difference. State the cost alongside any
result: `HNRNPC-Hela-iCLIP` was the only **assay-matched** column (THRAP3 is itself iCLIP),
and HNRNPC agreeing with *itself* across both cell line **and** assay was the strongest
finding of the earlier runs. HNRNPC keeps `HepG2-HNRNPC-eCLIP` and `K562-HNRNPC-eCLIP`, so it
is still an internal control — but a within-assay one, and the cross-assay argument no longer
applies. `ELAVL1-PARCLIP` and `TRA2B-MDA231-iCLIP` are excluded too. Point `-s` back at the
top-level `RBPeekSamplesheet.tsv` to restore all of them.

## Caveats when reading the results

- **Cell line is confounded with RBP.** THRAP3 is HEK293; the panel is HepG2/K562 eCLIP,
  HeLa iCLIP and PAR-CLIP. There is no HEK293 panel member. Rank RBPs against each other
  rather than reading absolute enrichment.
- **The peak near -5 is half of a symmetric +/-5 doublet, and it is an artefact.**
  `compute_counts_for_protein` attributes a panel peak's *entire* score to a single offset:
  its genomic **start** (`intersect_inference_bed.py:345,356`). That is exactly right for the
  1 nt crosslink sites the script was designed around, but our panel is Clippy *peaks*
  (mean width ~10 bp), so a peak centred on a THRAP3 anchor lands at offset `-w/2`. The
  minus-strand flip at `:357` then sends minus-strand loci to `+w/2`, because a peak's
  genomic start is its 3' end on that strand. Measured on 89,240 real intersections:
  **+ strand median -5.0, - strand median +5.0**, near-equal counts. The metaprofile
  x-axis is therefore a distribution of *peak-start positions*, not of crosslinks, and
  nothing at nucleotide resolution should be read from it.
  **Fixed by `--panel-anchor midpoint`**, now set in `run_thrap3_intersect.sbatch`. Scoring
  each panel interval at its midpoint collapses both modes onto 0: measured strand
  separation drops from 10.0 nt to 0.0, and the fraction of loci peaking within +/-2 nt
  goes from 7% to 31%. The flag is a no-op for 1 nt crosslink input, so it is safe to
  leave on if the panel is later switched to `../CLIP/*-merged-xls/*.xl.bed.gz` (268 files),
  which remains the better long-term input for nucleotide-resolution questions.

- **Read positional structure off `metaprofile.png`.** It now draws its top 10 from the same
  enrichment ranking that picks the heatmap's 20 columns, so the two figures agree on which
  proteins matter, and the separate `protein_nt_metaprofile_heatmap.png` it used to need has
  been removed. Do **not** read positional structure off the `-i/--inspect-protein` plot: it
  is one row per *locus*, hierarchically clustered with average linkage on sparse data, which
  chains — a 2-way cut of a 4,000-row reproduction split 3,999 vs 1, so its apparent row
  groups are not real sub-populations. These runs no longer pass `-i`.
- **The heatmap row filter is now a percentile**, `--heatmap-min-support-percentile 10`, not
  an absolute cut. Row support is the summed score of every overlapping panel interval across
  every column, so it has no intrinsic scale: median row support was 1,594 on the exonic
  subset against 296 on the intronic one, and the flat 200 used previously therefore dropped
  11.3% of exonic loci against 43.5% of intronic ones. The percentile also always drops
  zero-support loci — 18.6% of the intronic set — so a run where the percentile resolves to 0
  cannot produce blank heatmap rows.
- **These runs pass `--no-clustering`**, so there is no `binf_heatmap_clusters.tsv` and no
  `metaprofile_cluster_C*.png`. `protein_ranking.tsv` is the replacement: every panel sample
  ranked by mean peak support with its `frac centred` score alongside, which is the table to
  read when weighing depth against central enrichment.
- **Clippy parameterisation differs across the panel.** THRAP3 and the PAR-CLIP set use
  `minHeightAdjust1.0_minPromAdjust1.0`; the eCLIP set uses `stdev1.0`. This affects the
  panel columns' peak-calling sensitivity relative to each other.
