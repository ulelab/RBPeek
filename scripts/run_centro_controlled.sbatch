#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_centro_ctrl
#SBATCH --time=1-00:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-centro-%x-%j.out

# Centrosome apex CLIP, with the APEX controls applied.
#
#   sbatch --job-name=centro_specific scripts/run_centro_controlled.sbatch specific
#   sbatch --job-name=centro_ntctrl   scripts/run_centro_controlled.sbatch ntcontrol
#
# Run BOTH. Neither is interpretable alone: the point is to compare how often each panel
# protein binds centrosome-SPECIFIC loci versus how often it binds the untargeted-APEX
# control loci. A protein that binds both equally is reporting APEX background, not
# centrosome biology.
#
# WHY THIS EXISTS
# The uncontrolled run (results/centro) was dominated by background. Of its 8,218 anchors,
# 69.4% also appear in a control replicate, and those control-shared loci are precisely the
# promiscuous ones: median 181 of 302 panel proteins bound versus 114 for the specific set,
# and 21.8% of them are 3/3-reproducible versus 1.2%. Subtracting controls leaves 2,106
# anchors, of which only THREE are 3/3-reproducible - i.e. essentially all of the highly
# reproducible centrosome signal is shared with the controls.
#
# READ THE OUTPUT ACCORDINGLY
# --protein-select enrichment ranks by positional centring, and that metric was already
# non-discriminating on the full 8,218-locus set: the top-20 spanned only 0.043 (THRAP3
# spanned 0.096), ranks 2-20 sat within ~2 standard errors of each other, and the same
# protein landed at ranks 1/179/227/254 across its four datasets. On 2,106 loci it will be
# noisier still. Do not read the top-20 here as a hit list.
#
# The real output is the DIFFERENTIAL between the two runs, computed from the two
# binf_summary.tsv files: for each protein, the fraction of loci it binds in the specific
# run minus the fraction it binds in the control run. That does not depend on positional
# precision, which is the part proximity labelling cannot deliver - biotin labels within a
# radius, so peak position reflects distance from the centrosome rather than any RBP's own
# footprint.

set -euo pipefail

SET="${1:-}"
case "$SET" in
  specific)  BINF=Centrosome/Centro_apex_specific_merged_min2rep_anchors.bed; OUT=results/centro_specific ;;
  ntcontrol) BINF=Centrosome/NT_apex_merged_min2rep_anchors.bed;              OUT=results/centro_ntcontrol ;;
  *) echo "usage: sbatch $0 <specific|ntcontrol>" >&2; exit 2 ;;
esac

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/build_inference_bed_from_peaks.py first." >&2
  exit 1
fi
echo "set: $SET"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# --protein-min-loci 100, the SAME absolute value for both runs. Matching it matters more
# than tuning it: the differential is only meaningful across proteins that cleared the gate
# in both runs, and the two sets differ in size (2,106 vs 4,284). 100 is low enough for the
# smaller set to retain most of the panel. The enrichment fraction is correspondingly
# imprecise (SE ~ +/-0.04 at n=100), which is another reason to read the differential rather
# than the rank.
# --protein-nt-heatmap was removed from the tool: the metaprofile now takes its top proteins
# from the same enrichment ranking that picks the heatmap's columns, so it already shows the
# heatmap's proteins as profiles. Clustering is KEPT here - the Centrosome analysis reads
# cluster membership downstream via plot_offset_distribution.py and
# plot_cluster_metaprofile_at_loci.py, both of which need binf_heatmap_clusters.tsv.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s RBPeekSamplesheet.tsv \
  --skip-merge \
  --table \
  -o "$OUT" \
  --window 100 \
  --nt-heatmap-window 50 \
  --panel-anchor midpoint \
  --protein-select enrichment \
  --enrichment-window 5 \
  --centrality-min-total 100 \
  --protein-min-loci 100 \
  --heatmap-min-support 200 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 10 \
  --tsne
