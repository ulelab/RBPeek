#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_centro
#SBATCH --time=1-00:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-centro-%j.out

# Centrosome apex CLIP vs the RBP panel.
#
#   sbatch scripts/run_centro_intersect.sbatch
#
# Source: Flow execution 110200382017739198 (Centrosome_CLIP), three Clippy peak
# replicates Centro_apex_R1..R3, downloaded to Centrosome/raw/.
#
# Assembly confirmed from the data rather than assumed: autosomes run to 22, the KI270442.1
# scaffold is human, and the largest chr1 coordinate (247.6Mb) fits GRCh38 (248.9Mb) and not
# mm39 (195.1Mb). Safe against the human panel.
#
# Inference BED built by scripts/build_inference_bed_from_peaks.py: strand-aware merge,
# >=2/3 replicate support, collapsed to 1 nt midpoint anchors -> 8,218 loci.
#
# Replicate depth is very uneven: R1 25,709 peaks, R2 12,414, R3 1,629 - a 16x spread.
# R3 is shallow but precise (83% of its peaks are recovered by R1 and by R2), so it is
# underpowered rather than noisy. The kept set is nonetheless dominated by R1 n R2
# (6,778 of 8,218); R3 mostly contributes the 1,276-locus 3/3 core. Treat the 3/3 tier
# (score column == 3) as the high-confidence subset.
#
# NOTE ON THE BAIT: "Centro_apex" is a subcellular compartment capture, not a single RBP
# pulldown. The question this run answers is therefore "which RBPs bind at centrosome-
# localised RNA sites", not "which RBPs partner with protein X". There is no self-column to
# exclude from the panel and no assay-matched control to nominate.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

BINF=Centrosome/Centro_apex_merged_min2rep_anchors.bed
OUT=results/centro

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/build_inference_bed_from_peaks.py first." >&2
  exit 1
fi
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# Settings match the THRAP3 runs - see scripts/run_thrap3_region.sbatch for the reasoning.
#
# The per-sample eligibility gates these runs were originally tuned with are gone: every
# sample is now ranked, and protein_ranking.tsv reports loci_with_signal, so a sample topping
# the ranking on a handful of loci is visible rather than silently filtered. Check that
# column before trusting anything near the top.
#
# -n 6 is kept HERE and not in the THRAP3 runs: the Centrosome analysis reads cluster
# membership downstream via plot_offset_distribution.py --clusters and
# plot_cluster_metaprofile_at_loci.py --clusters, both of which need binf_heatmap_clusters.tsv.
#
# No exon/intron split for this analysis, by request - combined regions only.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s RBPeekSamplesheet.tsv \
  -o "$OUT" \
  --window 100 \
  --panel-anchor midpoint \
  --protein-select enrichment \
  --enrichment-window 5 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --top-proteins 20 \
  --support-pct 10 \
  --n-clusters 6 \
  --tsne \
