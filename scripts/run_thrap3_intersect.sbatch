#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3
#SBATCH --time=1-00:00:00
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%j.out

# THRAP3 vs the eCLIP panel, all anchors (no region split).
#
#   sbatch scripts/run_thrap3_intersect.sbatch [enrichment|total]
#
# THRAP3 is the INFERENCE BED (rows), not a panel column: it is absent from the samplesheet,
# so there is no self-comparison to exclude. Each row is a reproducible THRAP3 binding site;
# each column is one eCLIP peak set. The question is which RBPs co-occupy THRAP3 sites.
#
# Inference BED built by scripts/build_thrap3_inference_bed.py from the four Flow Clippy peak
# calls (project HA_THRAP3_CLIP 788995297969977723, CLIP-Seq v1.7, GRCh38): strand-aware
# merge, >=2/4 replicate support, collapsed to 1 nt midpoint anchors.
#
# Carried into interpretation: THRAP3 is HEK293 iCLIP, the panel is HepG2/K562 eCLIP. Cell
# line is confounded with RBP, so rank RBPs against each other rather than reading absolute
# enrichment. HNRNPC appears as both HepG2 and K562 eCLIP, so its two columns are an internal
# consistency check - within one assay, since the eCLIP sheet has no iCLIP.

set -euo pipefail

SELECT="${1:-enrichment}"
if [[ "$SELECT" != "enrichment" && "$SELECT" != "total" ]]; then
  echo "usage: sbatch $0 [enrichment|total]" >&2
  exit 2
fi

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

BINF=THRAP3/THRAP3_merged_min2rep_anchors.bed
OUT=results/thrap3
[[ "$SELECT" == "total" ]] && OUT="${OUT}_total"

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/build_thrap3_inference_bed.py first." >&2
  exit 1
fi
echo "ranking: $SELECT   -> $OUT"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# See scripts/run_thrap3_region.sbatch for why --panel-anchor midpoint is required here, what
# the eCLIP-only samplesheet costs, and what --support-pct 10 does. Add -n 6 for k-means
# clusters, binf_heatmap_clusters.tsv and per-cluster metaprofiles.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  -o "$OUT" \
  --window 100 \
  --panel-anchor midpoint \
  --protein-select "$SELECT" \
  --enrichment-window 5 \
  --top-proteins 20 \
  --support-pct 10 \
  --tsne
