#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3_region
#SBATCH --time=1-00:00:00
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%x-%j.out

# THRAP3 co-binding, split by transcript region.
#
#   sbatch --job-name=thrap3_exonic  scripts/run_thrap3_region.sbatch exonic
#   sbatch --job-name=thrap3_intronic scripts/run_thrap3_region.sbatch intronic
#
# (--job-name is worth passing: it lands in the %x of the output filename, so the two logs
# are told apart at a glance rather than by job id.)
#
# Same inference loci, same panel, same settings as run 52008929 - the ONLY difference
# between the two jobs is which subset of the THRAP3 anchors is used. That is what makes
# the two enrichment rankings comparable: any divergence is region, not parameters.
#
# Build the subsets once before submitting either job:
#   python3 scripts/split_inference_bed_by_region.py \
#     -b THRAP3/THRAP3_merged_min2rep_anchors.bed -g "$GTF" -o THRAP3
#
# Memory is raised to 90GB: the panel now includes four crosslink-level BEDs, one of which
# is ~550MB, far denser than the Clippy peak files the 60GB runs used.

set -euo pipefail

REGION="${1:-}"
if [[ "$REGION" != "exonic" && "$REGION" != "intronic" ]]; then
  echo "usage: sbatch $0 <exonic|intronic>" >&2
  exit 2
fi

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

BINF="THRAP3/THRAP3_merged_min2rep_anchors_${REGION}.bed"
OUT="results/thrap3_${REGION}"

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/split_inference_bed_by_region.py first." >&2
  exit 1
fi
echo "region: $REGION"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# Settings carried over unchanged from run 52008929 so the two regions stay comparable, and
# so both stay comparable to the whole-locus run. See scripts/run_thrap3_intersect.sbatch
# for why each one is set the way it is.
#
# --heatmap-min-support 200 is carried over, but WATCH IT. Row support sums the score column
# across every panel entry, and the panel just gained four crosslink-level BEDs that are far
# denser than the Clippy peak files (one is ~550MB against ~1MB typical). Row sums will rise,
# so 200 will be more permissive than the 78.3% it kept on run 52008929. The log prints the
# kept percentage - if it approaches 100%, raise this before reading the heatmap.
#
# --nt-heatmap-window 50 is new and is DISPLAY ONLY: --window stays at 100, so every
# statistic, the metaprofile, the enrichment ranking and the summary table are identical to
# what they would be uncropped. It only zooms both nucleotide heatmaps onto the centre.
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
  --protein-nt-heatmap \
  --protein-select enrichment \
  --enrichment-window 5 \
  --centrality-min-total 100 \
  --protein-min-loci 500 \
  --heatmap-min-support 200 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 10 \
  --tsne \
  -i HNRNPC-Hela-iCLIP
