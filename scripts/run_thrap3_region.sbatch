#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_thrap3_region
#SBATCH --time=1-00:00:00
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-thrap3-%x-%j.out

# THRAP3 co-binding analysis for one transcript region.
#
#   sbatch --job-name=thrap3_exonic   scripts/run_thrap3_region.sbatch exonic
#   sbatch --job-name=thrap3_intronic scripts/run_thrap3_region.sbatch intronic
#
# The locus sets and normalisation regions are generated beforehand with
#
#   python3 scripts/split_inference_bed_by_region.py \
#     -b THRAP3/THRAP3_merged_min2rep_anchors.bed -g <annotation.gtf> -o THRAP3 --drop-chrM
#
# Site-specific paths can be overridden through the environment:
#   PANEL_DIR  directory against which the samplesheet paths are resolved
#   GENOME     chromosome sizes file
#   CONDA_SH   conda initialisation script

set -euo pipefail

REGION="${1:-}"
if [[ "$REGION" != "exonic" && "$REGION" != "intronic" ]]; then
  echo "usage: sbatch $0 <exonic|intronic>" >&2
  exit 2
fi

PANEL_DIR="${PANEL_DIR:-../CLIP}"
GENOME="${GENOME:-/camp/home/jonesm6/home/shared/genomes/hg38/hg38.genome}"
CONDA_SH="${CONDA_SH:-/camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh}"

cd "${SLURM_SUBMIT_DIR:-.}"

source "$CONDA_SH"
conda activate rbpeek

BINF="THRAP3/THRAP3_merged_min2rep_anchors_${REGION}.bed"
NORM="THRAP3/regions_${REGION}.bed"
OUT="results/thrap3_${REGION}"

for f in "$BINF" "$NORM" "$GENOME"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: $f not found" >&2
    exit 1
  fi
done
echo "region: $REGION   -> $OUT"
echo "inference loci: $(wc -l < "$BINF")"

mkdir -p "$OUT"

python3 scripts/intersect_inference_bed.py \
  -x "$PANEL_DIR" \
  -b "$BINF" \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  --norm-bed "$NORM" \
  --genome "$GENOME" \
  -o "$OUT" \
  --window 100 \
  --central-window 10 \
  --support-pct 20 \
  --tsne
