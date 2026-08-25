# Centrosome apex CLIP vs the RBP panel

Which RBPs bind at centrosome-localised RNA sites, and with what profile shape.

## Source

Flow execution [110200382017739198](https://app.flow.bio/executions/110200382017739198)
(Centrosome_CLIP). Three Clippy peak replicates, `Centro_apex_R1..R3`, in `raw/`.

**Assembly verified from the data, not assumed** — the Flow API token had expired, so
provenance was checked directly: autosomes run to 22 (mouse stops at 19), the `KI270442.1`
scaffold is human, and the largest chr1 coordinate (247.6 Mb) fits GRCh38 (248.9 Mb) rather
than mm39 (195.1 Mb). Human, safe against the human panel.

**The bait is a compartment, not a protein.** "Centro_apex" is a subcellular capture, so
there is no self-column to exclude from the panel and no assay-matched control to nominate
as `--inspect-protein`. Read results as "RBPs at centrosome-localised sites", not "RBPs
partnering with protein X".

## Replicates

| replicate | peaks | notes |
|---|---:|---|
| R1 | 25,709 | deepest |
| R2 | 12,414 | |
| R3 | 1,629 | 16x shallower than R1 |

Pairwise recovery (% of row's peaks overlapped by column, strand-aware):

| | ->R1 | ->R2 | ->R3 |
|---|---|---|---|
| **R1** | – | 31.5% | 5.4% |
| **R2** | 65.0% | – | 11.0% |
| **R3** | **83.7%** | **82.7%** | – |

R3 is shallow but **precise** — 83% of its peaks are recovered by each of the others, so it
is underpowered rather than noisy. A clean depth ladder, each replicate nested in the deeper
ones.

## Inference BED

```bash
python3 scripts/build_inference_bed_from_peaks.py \
  -p Centrosome/raw/Centro_apex_R1_*_Peaks.bed \
     Centrosome/raw/Centro_apex_R2_*_Peaks.bed \
     Centrosome/raw/Centro_apex_R3_*_Peaks.bed \
  -l Centro_apex --min-reps 2 -o Centrosome
```

Strand-aware merge, `>=2/3` replicate support, collapsed to 1 nt midpoint anchors.
**8,218 loci**, max merged width 59 bp so no over-chaining.

| support | regions | |
|---|---:|---|
| 1/3 | 21,675 | 72.5% |
| 2/3 | 6,942 | 23.2% |
| 3/3 | 1,276 | 4.3% |

Composition of the kept set: `R1,R2` 6,778 · `R1,R2,R3` 1,276 · `R1,R3` 90 · `R2,R3` 74.

**The kept set is essentially R1 n R2.** R3's distinct contribution is the 1,276-locus 3/3
core. The score column carries replicate support, so `score == 3` selects that
high-confidence subset for any follow-up.

## Running

```bash
sbatch scripts/run_centro_intersect.sbatch
```

No exon/intron split for this analysis — combined regions only.

Settings match THRAP3 run 52008929 except `--protein-min-loci 300` (not 500). That gate is
an absolute count of loci where a protein has signal, and this set is 8,218 loci against
THRAP3's 29,018 — 500 would demand signal at 6.1% of loci here versus 1.7% there, excluding
much more of the panel. 300 asks for 3.7%; the standard error on an enrichment fraction near
0.30 widens only from +/-0.021 to +/-0.027, well inside the ~0.10 spread the top 20 covers.

## What to check in the log

- **`Enrichment ranking: N of 302 proteins excluded`** — if N is large, `--protein-min-loci`
  is biting and should come down further.
- **row filter kept %** — `--heatmap-min-support 200` sums support across all 302 columns,
  so it does not scale with locus count. If it keeps ~everything, raise it.
- **`-i HNRNPC-Hela-iCLIP`** is carried over from the THRAP3 runs for continuity only; there
  is no bait-matched reason for it here. Swap it once the ranking shows what is actually
  enriched.
