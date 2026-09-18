#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3_region
#SBATCH --time=1-00:00:00
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%x-%j.out

# THRAP3 co-binding, split by transcript region.
#
#   sbatch --job-name=thrap3_exonic   scripts/run_thrap3_region.sbatch exonic
#   sbatch --job-name=thrap3_intronic scripts/run_thrap3_region.sbatch intronic
#
# Build the locus subsets and the normalisation regions once before submitting 
#
#   python3 scripts/split_inference_bed_by_region.py \
#     -b THRAP3/THRAP3_merged_min2rep_anchors.bed -o THRAP3 --drop-chrM


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
python -c "import matplotlib, seaborn, sklearn; print('ok')"

BINF="THRAP3/THRAP3_merged_min2rep_anchors_${REGION}.bed"
NORM="THRAP3/regions_${REGION}.bed"
OUT="results/thrap3_${REGION}"

for f in "$BINF" "$NORM"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: $f missing. Run scripts/split_inference_bed_by_region.py first (see header)." >&2
    exit 1
  fi
done
echo "region: $REGION   -> $OUT"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  --norm-bed "$NORM" \
  -o "$OUT" \
  --window 100 \
  --support-pct 20 \
  --save-matrix 10 \
  --tsne
