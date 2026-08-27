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
# Build the subsets once before submitting either job (the GTF path is now the default,
# so -g is only needed to override it):
#   python3 scripts/split_inference_bed_by_region.py \
#     -b THRAP3/THRAP3_merged_min2rep_anchors.bed -o THRAP3
#
# Default annotation: gencodev39_annotation/filtered.gencode.v39.main.annotation.gtf, the
# same flavour Clippy was run against, and chr-prefixed like the inference BED.
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
# --protein-min-loci-frac 0.0248 REPLACES the flat --protein-min-loci 500. The gate is
# load-bearing either way - the enrichment fraction is trivially 1.0 for a protein whose
# signal touches a single locus - but an ABSOLUTE gate cannot be shared by two runs of
# different size. Exonic has 20,158 loci and intronic 8,666, so a flat 500 demanded signal at
# 2.5% of loci in one and 5.8% in the other. Runs 52359304/52359306 show the damage: 38 of
# 302 proteins excluded from exonic against 90 from intronic, so a protein binding ~4% of loci
# was ranked in one run and invisible in the other purely on set size - and in the output,
# exclusion is indistinguishable from absence.
#
# 0.0248 is 500/20158: what 500 meant on the exonic set. It resolves to 500 for exonic and 215
# for intronic, gating both identically in relative terms. The cost is precision on the
# smaller set - the standard error on an enrichment fraction near 0.30 widens from +/-0.021 at
# n=500 to +/-0.031 at n=215. The run log prints the resolved value.
#
# BOTH REGIONS MUST BE RESUBMITTED regardless: runs 52359304/52359306 predate the chromosome
# harmonisation fix, so K562-ELAVL1-rep1/rep2 and TRA2B-MDA231-rep1..rep3 contributed nothing
# to either. Their logs have no "Chromosome naming:" line; the new ones should.
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
  --protein-min-loci-frac 0.0248 \
  --heatmap-min-support 200 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --cluster-metaprofiles \
  --n-clusters 6 \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 10 \
  --tsne \
  -i HNRNPC-Hela-iCLIP
