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

# Settings carried over from THRAP3 run 52008929 - see scripts/run_thrap3_intersect.sbatch
# for the reasoning behind each. Two deliberate differences:
#
# --protein-min-loci 300 rather than 500. The gate is an ABSOLUTE count of loci where a
# protein has signal, but this inference set is 8,218 loci against THRAP3's 29,018, so 500
# would demand signal at 6.1% of loci here versus 1.7% there, excluding far more of the
# panel. 300 asks for 3.7%. The cost is precision on the enrichment fraction: the standard
# error on a fraction near 0.30 goes from +/-0.021 at n=500 to +/-0.027 at n=300, still well
# inside the ~0.10 spread the top 20 usually covers. Raise it back to 500 if the selected
# columns look thin.
#
# No exon/intron split for this analysis, by request - combined regions only.
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
  --protein-min-loci 300 \
  --heatmap-min-support 200 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 10 \
  --tsne \
  -i HNRNPC-Hela-iCLIP
