#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_intersect
#SBATCH --time=1-00:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=20

set -euo pipefail

python3 scripts/intersect_inference_bed.py \
  --xldir ../CLIP \
  --bed decoys/Decoys_proteincoding_splicescores.bed \
  --samplesheet RBPeekSamplesheet.tsv \
  --skip-merge \
  --table \
  --outdir results
