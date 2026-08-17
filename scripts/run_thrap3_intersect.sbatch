#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3
#SBATCH --time=1-00:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%j.out

# THRAP3 vs the 297-sample RBP panel.
#
# THRAP3 is the INFERENCE BED (rows), not a panel column: it is absent from
# RBPeekSamplesheet.tsv, so there is no self-comparison to exclude. Each row is a
# reproducible THRAP3 binding site; each column is one of the 297 Clippy peak sets
# under ../CLIP. The question the run answers is "which RBPs co-occupy THRAP3 sites,
# and with what profile shape".
#
# Inference BED built by scripts/build_thrap3_inference_bed.py from the four Flow
# Clippy peak calls (project HA_THRAP3_CLIP 788995297969977723, CLIP-Seq v1.7, GRCh38):
# strand-aware merge, >=2/4 replicate support, collapsed to 1 nt midpoint anchors.
#
# Caveat carried into interpretation: THRAP3 is HEK293 iCLIP, while the panel is
# HepG2/K562 eCLIP, HeLa iCLIP and PAR-CLIP. Cross-cell-line differences are
# confounded with cross-RBP differences, so rank RBPs against each other rather than
# reading absolute enrichment. HNRNPC-Hela-iCLIP is the assay-matched control and is
# used as --inspect-protein for that reason; HNRNPC also appears in four independent
# assay/cell combinations, so its four columns act as an internal consistency check.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

BINF=THRAP3/THRAP3_merged_min2rep_anchors.bed
OUT=results/thrap3

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/build_thrap3_inference_bed.py first." >&2
  exit 1
fi
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# --panel-anchor midpoint is REQUIRED for this panel. RBPeekSamplesheet.tsv points at
# Clippy *peaks* (mean width ~11 bp), and a panel interval's whole score lands on a single
# offset. Scoring at the interval start therefore puts a centred peak at -w/2 on + strand
# loci and +w/2 on - strand loci, splitting real central signal into a +/-5 doublet with a
# hole at 0. Measured on this panel: separation 10.0 nt with start, 0.0 nt with midpoint.
# If the panel is ever switched to the 1 nt crosslink files (../CLIP/*-merged-xls/), the
# flag becomes a no-op rather than a problem, since start == midpoint at width 1.
#
# --protein-select centrality replaces depth ranking. Ranking by summed total_overlaps put
# peak-rich datasets at the top whether or not their binding relates to THRAP3: in a
# controlled test, a decoy displaced 40 nt and a decoy with jittered positions both ranked
# above genuine replicates purely on depth. Under centrality every genuine replicate
# outranked both decoys, including one at a quarter of their depth, and the 40 nt decoy
# scored r=-0.003.
#
# --heatmap-min-support 200, deliberately low. This filter is NOT a peak count: row support
# is the summed SCORE column of every overlapping panel interval across all 297 columns and
# the whole +/-100nt window, so its scale is arbitrary (median 1,053 on the last run).
#
# It also correlates with THRAP3 reproducibility - median row support is 821 at 2/4 reps and
# 2,324 at 4/4 - so a high threshold silently culls the least reproducible loci. Retention
# of 2/4-rep vs 4/4-rep loci: 87%/87% at 50, 76%/82% at 200, 61%/74% at 500, 29%/53% at
# 2000. 200 keeps the bias mild while still dropping near-empty rows from the k-means.
#
# Raising it does NOT fix the dark heatmap. That comes from logistic_scale() centring on the
# matrix median, which is ~0 for a sparse matrix, so every empty cell renders at exactly 0.5
# and the whole palette below 0.5 goes unused - which is why the last run's colourbar
# started at 0.5. --heatmap-scale percentile is the actual fix: it log1p's the matrix and
# scales against the 99th percentile of non-zero cells, so empty cells sit at 0 and the full
# palette carries signal. On a sparse test matrix (65% empty): logistic gave range
# [0.50,1.00] with empty cells at 0.50; percentile gave [0.00,1.00] with empty cells at 0.00
# and the median non-zero cell at 0.53. The log1p step matters - support counts are
# heavy-tailed, and without it the median non-zero cell rendered at 0.11, darker than before.
#
# Nothing is lost either way: binf_summary.tsv and binf_heatmap_clusters.tsv keep all 29,018
# loci with a passes_heatmap_filter flag, and the metaprofile, the proteins x nt heatmap and
# the centrality ranking are all computed over every locus, unfiltered.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s RBPeekSamplesheet.tsv \
  --skip-merge \
  --table \
  -o "$OUT" \
  --window 100 \
  --panel-anchor midpoint \
  --protein-nt-heatmap \
  --protein-select centrality \
  --centrality-sigma 5 \
  --centrality-min-total 100 \
  --heatmap-min-support 200 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 10 \
  --tsne \
  -i HNRNPC-Hela-iCLIP
