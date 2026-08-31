#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3_region
#SBATCH --time=1-00:00:00
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%x-%j.out

# THRAP3 co-binding, split by transcript region.
#
#   sbatch --job-name=thrap3_exonic            scripts/run_thrap3_region.sbatch exonic
#   sbatch --job-name=thrap3_intronic          scripts/run_thrap3_region.sbatch intronic
#   sbatch --job-name=thrap3_intronic_total    scripts/run_thrap3_region.sbatch intronic total
#
# Second argument picks the ranking: enrichment (default) or total. It only changes WHICH
# panel samples reach the figures - protein_ranking.tsv carries both rankings either way, so
# two runs differing only in this argument are directly comparable row for row.
#
# --job-name is worth passing: it lands in the %x of the output filename, so the logs are
# told apart at a glance rather than by job id.
#
# Build the region subsets once before submitting (the GTF path is the script's default):
#   python3 scripts/split_inference_bed_by_region.py \
#     -b THRAP3/THRAP3_merged_min2rep_anchors.bed -o THRAP3
#
# Memory is 90GB: the panel includes crosslink-level BEDs far denser than the Clippy peak
# files the older 60GB runs used.

set -euo pipefail

REGION="${1:-}"
SELECT="${2:-enrichment}"
if [[ "$REGION" != "exonic" && "$REGION" != "intronic" ]]; then
  echo "usage: sbatch $0 <exonic|intronic> [enrichment|total]" >&2
  exit 2
fi
if [[ "$SELECT" != "enrichment" && "$SELECT" != "total" ]]; then
  echo "usage: sbatch $0 <exonic|intronic> [enrichment|total]" >&2
  exit 2
fi

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

BINF="THRAP3/THRAP3_merged_min2rep_anchors_${REGION}.bed"
OUT="results/thrap3_${REGION}"
[[ "$SELECT" == "total" ]] && OUT="${OUT}_total"

if [ ! -s "$BINF" ]; then
  echo "ERROR: $BINF missing. Run scripts/split_inference_bed_by_region.py first." >&2
  exit 1
fi
echo "region: $REGION   ranking: $SELECT   -> $OUT"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

# --panel-anchor midpoint is REQUIRED for this panel. The samplesheet points at Clippy PEAKS
# (mean width ~11 bp), and a panel interval's whole score lands on one offset. Scoring at the
# interval start puts a centred peak at -w/2 on + strand loci and +w/2 on - strand loci,
# splitting real central signal into a +/-5 doublet with a hole at 0. Measured on this panel:
# 10.0 nt separation with start, 0.0 with midpoint.
#
# -s THRAP3/RBPeekSamplesheet_eCLIP.tsv is the 224 eCLIP columns of 299. Assay is uniform, so
# a difference between two columns cannot be an assay difference. The cost: HNRNPC-Hela-iCLIP
# was the only assay-matched column (THRAP3 is itself iCLIP) and HNRNPC agreeing with itself
# across cell line AND assay was the strongest earlier result. HNRNPC keeps HepG2 and K562, so
# it is still a control, but a within-assay one. Point -s at RBPeekSamplesheet.tsv to restore.
#
# --support-pct 10 drops the bottom decile of loci by summed support, and always drops loci
# with zero support. An absolute threshold cannot be shared by these two runs: median row
# support is 1,594 exonic against 296 intronic, so the flat 200 used previously dropped 11.3%
# of exonic loci against 43.5% of intronic ones. Expected here: exonic threshold 166 / 10.0%
# dropped; intronic threshold 0, zero clause fires, 18.6% dropped.
#
# No -n, so loci are ordered by total support and no clusters are written. Add -n 6 to get
# k-means groups, binf_heatmap_clusters.tsv and the per-cluster metaprofiles back.
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
