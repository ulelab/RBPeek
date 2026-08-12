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
# HepG2/K562 eCLIP, HeLa iCLIP and PAR-CLIP. Cross-cell-line differences are
# confounded with cross-RBP differences, so rank RBPs against each other rather than
# reading absolute enrichment. HNRNPC-Hela-iCLIP is the assay-matched control and is
# used as --inspect-protein for that reason; HNRNPC also appears in four independent
# assay/cell combinations, so its four columns act as an internal consistency check.

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

python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s RBPeekSamplesheet.tsv \
  --skip-merge \
  --table \
  -o "$OUT" \
  --window 100 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 60 \
  --metaprofile-top-proteins 10 \
  --tsne \
  -i HNRNPC-Hela-iCLIP
