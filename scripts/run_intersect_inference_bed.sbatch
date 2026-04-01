#!/usr/bin/env bash
#SBATCH --job-name=rbpeek_intersect
#SBATCH --time=1-00:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=20

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate rbpeek

which python
python -c "import matplotlib, seaborn, pybedtools, sklearn; print('ok')"

# Default row clustering is DBSCAN (--cluster-method dbscan). For average linkage + fixed k clusters:
#   add: --cluster-method hierarchical --n-clusters 4
# Optional: --cluster-top-proteins 30 (default), --dbscan-eps / --dbscan-min-samples to tune DBSCAN

python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b decoys/Decoys_proteincoding_splicescores.bed \
  -s RBPeekSamplesheet.tsv \
  --skip-merge \
  --table \
  -o results \
  --cluster-method dbscan \
  --cluster-metaprofiles \
  --tsne
