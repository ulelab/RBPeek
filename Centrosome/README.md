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
to nominate as a bait-matched control. Read results as "RBPs at centrosome-localised sites", not "RBPs
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

Settings match the THRAP3 runs. NOTE: the per-sample eligibility gates these runs were
tuned with (`--protein-min-loci`, `--centrality-min-total`) no longer exist — every sample is
now ranked, and `protein_ranking.tsv` reports `loci_with_signal` so a top-ranked sample
resting on a handful of loci is visible rather than silently filtered. The old note read:
That gate was
an absolute count of loci where a protein has signal, and this set is 8,218 loci against
THRAP3's 29,018 — 500 would demand signal at 6.1% of loci here versus 1.7% there, excluding
much more of the panel. 300 asks for 3.7%; the standard error on an enrichment fraction near
0.30 widens only from +/-0.021 to +/-0.027, well inside the ~0.10 spread the top 20 covers.

## What to check in the log

- **`loci_with_signal` in `protein_ranking.tsv`** — check it for anything near the top of
  the ranking. The enrichment fraction is trivially 1.0 for a sample touching one locus.
- **row filter kept %** — `--support-pct 10` drops the bottom decile of loci by summed
  support, plus every zero-support locus. It replaced an absolute threshold, which could not
  be shared between runs over different locus sets.
- **`-i HNRNPC-Hela-iCLIP`** is carried over from the THRAP3 runs for continuity only; there
  is no bait-matched reason for it here. Swap it once the ranking shows what is actually
  enriched.

---

## Controls (added later) — and what they revealed

Two APEX controls, three replicates each, in `raw/`:

| control | peaks (R1/R2/R3) | what it controls for |
|---|---|---|
| `NT_apex` | 7,542 / 3,555 / 6,078 | APEX active but untargeted — labelling anywhere in the cell |
| `No_apex` | 2,923 / **29** / 8,109 | no enzyme — endogenous biotinylation and non-specific pulldown |

`No_apex_R2` has 29 peaks and is effectively a failed library; treat `No_apex` as R1+R3.

### The uncontrolled run was mostly background

Of the 8,218 anchors in `Centro_apex_merged_min2rep_anchors.bed`, **69.4% overlap a control
peak** (NT 61.6%, No_apex 50.5%). And the control-shared loci are exactly the promiscuous
ones:

| | n | median panel proteins bound | median total signal | % 3/3-rep |
|---|---:|---:|---:|---:|
| control-shared | 5,706 | **181** of 302 | 10,961 | **21.8%** |
| centrosome-specific | 2,512 | **114** | 3,933 | **1.2%** |

Subtracting controls at the merged-region level leaves **2,106 anchors**, of which only
**three** are 3/3-reproducible — against 1,276 before subtraction. Essentially all of the
highly reproducible centrosome signal is shared with the controls.

That is the headline: the promiscuity you spotted in the heatmap is APEX background labelling
abundant, accessible transcripts, not multivalent binding at centrosome-specific sites.

### Rebuild

```bash
R=_genome.clippy._rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed
python3 scripts/build_inference_bed_from_peaks.py \
  -p raw/Centro_apex_R{1,2,3}$R -s raw/NT_apex_R{1,2,3}$R raw/No_apex_R{1,2,3}$R \
  -l Centro_apex_specific --min-reps 2 -o Centrosome
python3 scripts/build_inference_bed_from_peaks.py \
  -p raw/NT_apex_R{1,2,3}$R -l NT_apex --min-reps 2 -o Centrosome
```

`--subtract` drops any merged region overlapping a control **before** the reproducibility
filter, so a region is judged non-specific on *any* control evidence. For a proximity bait
the background is the thing being controlled for, so that is the conservative direction.

### Running the controlled comparison

```bash
sbatch --job-name=centro_specific scripts/run_centro_controlled.sbatch specific
sbatch --job-name=centro_ntctrl   scripts/run_centro_controlled.sbatch ntcontrol
```

Run **both** — neither is interpretable alone. Then:

```bash
python3 scripts/compare_binding_frequency.py \
  -a results/centro_specific/binf_summary.tsv  --label-a centro_specific \
  -b results/centro_ntcontrol/binf_summary.tsv --label-b NT_control \
  -o results/centro_specific/binding_frequency_vs_control.tsv
```

### Why the differential, and not the enrichment rank

`--protein-select enrichment` ranks by positional centring, and that metric does not
discriminate on this assay. On the full 8,218-locus set the top-20 spanned only 0.043
(THRAP3 spanned 0.096), ranks 2–20 sat within ~2 standard errors of one another, and the
same protein landed at wildly different ranks across its own datasets:

| protein | ranks across its datasets |
|---|---|
| TIA1 | 1, 179, 227, 254 |
| HNRNPC | 7, 134, 144, 264 |
| HNRNPL | 4, 70, 145 |

On 2,106 loci it will be noisier still. **Do not read the top-20 as a hit list.**

What *does* survive is which loci a protein binds: after regressing out the per-locus mean,
same-protein pairs still correlate at **+0.264** against **−0.022** for different-protein
pairs. So compare binding *frequency* between the specific and control sets, which needs no
positional precision — the thing proximity labelling cannot give, since biotin labels within
a radius and peak position reflects distance from the compartment rather than an RBP
footprint.
