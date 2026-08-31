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
# HepG2/K562 eCLIP in the eCLIP-only sheet this now uses. Cross-cell-line differences are
# confounded with cross-RBP differences, so rank RBPs against each other rather than
# reading absolute enrichment. HNRNPC appears as both HepG2 and K562 eCLIP, so its two
# columns still act as an internal consistency check - but within one assay now, since
# HNRNPC-Hela-iCLIP is not in the eCLIP sheet.

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
# --protein-select enrichment, NOT centrality. Both are depth-free, but centrality scores
# Pearson r of the mean profile against a centred Gaussian, and once --panel-anchor midpoint
# centres everything that stops discriminating on position and starts discriminating on how
# tidy a dataset's baseline is. On run 52005272 that put PAR-CLIP at 2.2x its share of the
# panel and pushed HepG2-BCLAF1-eCLIP - THRAP3's known complex partner and the tallest curve
# in the metaprofile - out of the top 20 entirely.
#
# enrichment instead asks, per locus, whether the protein's strongest binding lands within
# --enrichment-window of the THRAP3 site, then takes the fraction of that protein's
# signal-bearing loci where it does. Recomputed on run 52005272's summary table this gives
# 19 eCLIP + 1 iCLIP (panel composition is 75/23/2 eCLIP/PAR-CLIP/iCLIP), with
# HepG2-HNRNPC-eCLIP and HNRNPC-Hela-iCLIP first and second and BCLAF1 back at 16th.
#
# --protein-min-loci 500 is load-bearing. The metric is trivially 1.0 for a protein whose
# signal touches one locus: ungated, the top of this ranking was WDR33-PARCLIP with a single
# locus and a summed score of 5.
#
# -s THRAP3/RBPeekSamplesheet_eCLIP.tsv restricts the panel to the 224 eCLIP columns of 299,
# dropping 68 PAR-CLIP and 7 iCLIP, so assay is uniform and a difference between columns
# cannot be an assay difference. State the cost alongside the results: HNRNPC-Hela-iCLIP was
# the only assay-matched column (THRAP3 is itself iCLIP) and HNRNPC agreeing with ITSELF
# across both cell line AND assay was the strongest result of the earlier runs. HNRNPC keeps
# HepG2 and K562, but that replication is now within-assay. TRA2B-MDA231-iCLIP goes too.
# Point -s back at RBPeekSamplesheet.tsv to restore both. This also makes the old
# --exclude-groups PARCLIP argument moot - PAR-CLIP is simply not in the sheet.
#
# --heatmap-min-support-percentile 10 replaces the absolute --heatmap-min-support 200. Row
# support is the summed SCORE of every overlapping panel interval across every column, so it
# has no intrinsic scale and an absolute threshold means different things in different runs -
# median row support was 1,594 on the exonic subset against 296 on the intronic one. The
# percentile also always drops zero-support loci, whatever it resolves to, so a locus set
# where more than 10% of loci are empty cannot end up with blank heatmap rows.
#
# Raising the filter never fixed a dark heatmap; --heatmap-scale percentile did. It log1p's
# the matrix and scales against the 99th percentile of non-zero cells, so empty cells sit at 0
# and the whole palette carries signal. logistic_scale() centred on the matrix median, which
# is ~0 when the matrix is sparse, so every empty cell rendered at exactly 0.5 and the palette
# below 0.5 went unused. On a sparse test matrix (65% empty): logistic gave [0.50,1.00] with
# empty cells at 0.50; percentile gives [0.00,1.00] with empty at 0.00 and the median non-zero
# cell at 0.53. The log1p step matters - support counts are heavy-tailed, and without it the
# median non-zero cell rendered at 0.11, darker than before.
#
# --no-clustering. k-means on binarised presence answers "what kinds of locus are there",
# a Centrosome question; here the question is which panel columns are enriched at THRAP3
# sites. It removes binf_heatmap_clusters.tsv and the per-cluster metaprofiles, orders the
# heatmap's locus axis by total support, and leaves the tSNE one colour. protein_ranking.tsv
# is the replacement: every panel sample ranked by mean peak support, carrying frac_centred
# alongside, so depth and central enrichment sit in one table instead of one being picked in
# advance.
#
# No -i and no --protein-nt-heatmap: the metaprofile now takes its top 10 from the SAME
# enrichment ranking that picks the heatmap's 20 columns, which was the only job the
# proteins x nt heatmap had.
#
# Nothing is lost by the row filter either way: binf_summary.tsv keeps all 29,018 loci, and
# the metaprofile, the ranking and the enrichment scores are computed over every locus,
# unfiltered.
# --metaprofile-top-proteins 20 matches --cluster-top-proteins 20, so the metaprofile and the
# heatmap now show the IDENTICAL set rather than the metaprofile showing the top half of it.
# At 10 they diverged in a way that read as a bug: enrichment rank and depth are only weakly
# related here (Spearman +0.35 over the 224 eligible columns), so the metaprofile's ranks 1-10
# are mostly shallow columns while the heatmap's brightest rows - BCLAF1, BUD13, SAFB - sit at
# enrichment ranks 16, 17 and 12 and were drawn on the heatmap but absent from the curves.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  --skip-merge \
  --table \
  -o "$OUT" \
  --window 100 \
  --panel-anchor midpoint \
  --protein-select enrichment \
  --enrichment-window 5 \
  --centrality-min-total 100 \
  --protein-min-loci 500 \
  --heatmap-min-support-percentile 10 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --no-clustering \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 20 \
  --tsne
