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
`HNRNPC-Hela-iCLIP` is used as `--inspect-protein` because it is the only **assay-matched**
column (iCLIP, like THRAP3).

## Caveats when reading the results

- **Cell line is confounded with RBP.** THRAP3 is HEK293; the panel is HepG2/K562 eCLIP,
  HeLa iCLIP and PAR-CLIP. There is no HEK293 panel member. Rank RBPs against each other
  rather than reading absolute enrichment.
- **Expect the metaprofile to centre near -5, not 0.** Panel signal is attributed at panel
  peak *start*, while anchors are region *midpoints*, and mean peak width is ~11 bp. Verified
  on a local smoke test (median `max_binding_offset` -4 to -5). It is a uniform convention
  shift across all 297 columns, not biology. `build_thrap3_inference_bed.py --anchor start`
  removes it if a 0-centred plot is preferred.
- **The heatmap row filter is hardcoded** at `sum >= 50` (`intersect_inference_bed.py:704`),
  not the 40 stated in the top-level README. It was tuned against the decoys BED; with 297
  columns over 29,018 THRAP3 loci most rows will pass, so it may need raising to keep the
  heatmap legible.
- **Clippy parameterisation differs across the panel.** THRAP3 and the PAR-CLIP set use
  `minHeightAdjust1.0_minPromAdjust1.0`; the eCLIP set uses `stdev1.0`. This affects the
  panel columns' peak-calling sensitivity relative to each other.
