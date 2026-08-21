#!/usr/bin/env bash
#SBATCH --job-name=clippy_new
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-clippy-%j.out

# Call Clippy peaks for the two genuinely new samples.
#
# Why this job exists: the RBPeek panel is Clippy PEAK calls, but the files added to
# /camp/lab/ulej/home/shared/RBPeek are crosslink (xl) BEDs. Peaks for these two samples
# do not exist anywhere under ../CLIP - checked against CLIPtoc.txt, which has zero
# ELAVL1 eCLIP and zero TRA2B/MDA231 entries of any kind. The only ELAVL1 peaks in the
# panel are PAR-CLIP (ELAVL1-PARCLIP), already present.
#
# Not called here, deliberately:
#   TRA2A_HepG2 / TRA2A_K562 - already in the panel as HepG2-TRA2A-eCLIP and
#   K562-TRA2A-eCLIP, called from the same underlying experiments.
#   The two _Mm files - mouse, and the inference BED is hg38.
#
# Parameters match local-nf-modules/clippy/main.nf and nextflow.config in the
# peak-benchmarking repo (-n 10 -x 1 -mg 5 -mb 5), which is also how the THRAP3 peaks
# themselves were called on Flow. Output filenames therefore carry the
# rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5 convention.
#
# NOTE the panel's eCLIP entries were called by an older Clippy whose filenames read
# stdev1.0 instead. That is a real parameterisation difference between these two new
# columns and the 223 eCLIP columns, affecting peak-calling sensitivity - worth
# remembering when comparing their enrichment ranks against the rest of the panel.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"

source /camp/home/jonesm6/home/users/jonesm6/software/miniconda3/etc/profile.d/conda.sh
conda activate peak_benchmarking   # the env that provides clippy
which clippy

# EDIT THESE TWO if the paths differ on your side.
GTF="${CLIPPY_GTF:-/camp/home/jonesm6/home/shared/genomes/hg38/Homo_sapiens_filtered.gtf}"
FAI="${CLIPPY_FAI:-/camp/home/jonesm6/home/shared/genomes/hg38/Homo_sapiens.GRCh38.fasta.fai}"

for f in "$GTF" "$FAI"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: not found: $f" >&2
    echo "Set CLIPPY_GTF / CLIPPY_FAI, or edit this script." >&2
    exit 1
  fi
done

OUTDIR=clippy_peaks
mkdir -p "$OUTDIR"

run_one () {
  local xl="$1" prefix="$2"
  if [ ! -s "$xl" ]; then
    echo "ERROR: missing input $xl" >&2
    exit 1
  fi
  echo "=== $prefix  <- $xl"
  clippy \
    -i "$xl" \
    -o "$OUTDIR/${prefix}" \
    -a "$GTF" \
    -g "$FAI" \
    -n 10 \
    -x 1 \
    -mg 5 \
    -mb 5 \
    -t "${SLURM_CPUS_PER_TASK:-12}" \
    -cf 16
}

run_one ELAVL1_K562_eCLIP_Hs.genome.xl.merged.bed   K562-ELAVL1-eCLIP
run_one TRA2B_MDA231_iCLIP_Hs.genome.xl.merged.bed  TRA2B-MDA231-iCLIP

echo
echo "=== peaks written ==="
ls -l "$OUTDIR"/*_Peaks.bed 2>/dev/null || ls -l "$OUTDIR"

echo
echo "Next: append these to RBPeekSamplesheet.tsv (paths resolve from -x ../CLIP),"
echo "substituting the real filenames printed above:"
for p in "$OUTDIR"/*_Peaks.bed; do
  [ -e "$p" ] || continue
  base=$(basename "$p")
  grp=$(echo "$base" | sed -E 's/_rollmean.*//')
  printf '  ../RBPeek/%s/%s\t%s\n' "$OUTDIR" "$base" "$grp"
done
