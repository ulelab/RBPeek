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
# --heatmap-min-support-percentile 10 REPLACES the absolute --heatmap-min-support 200, for
# exactly the reason --protein-min-loci-frac replaced --protein-min-loci below. Row support
# has no intrinsic scale: it is the summed SCORE of every overlapping panel interval across
# every column, so the same number means different things in two runs. Measured on runs
# 52359304/52359306, median row support was 1,594 exonic against 296 intronic - a 5.4x gap -
# so a flat 200 dropped 11.3% of exonic loci and 43.5% of intronic ones. The intronic heatmap
# was effectively a top-56% subset of its loci while the exonic one was a top-89% subset, and
# the two were being read side by side.
#
# The percentile filter also ALWAYS drops zero-support loci, whatever the percentile resolves
# to. That clause is load-bearing here: 18.6% of intronic loci have no support from any panel
# column, so the 10th percentile of the intronic set is literally 0 and a bare percentile
# would have kept ~750 blank heatmap rows. Expected: exonic threshold 166, 10.0% dropped;
# intronic threshold 0 with the zero clause firing, 18.6% dropped. The residual asymmetry is
# now real emptiness rather than an artefact of scale.
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
# -s THRAP3/RBPeekSamplesheet_eCLIP.tsv restricts the panel to the 224 eCLIP columns of 299,
# dropping 68 PAR-CLIP and 7 iCLIP. Assay is then uniform, so a difference between two columns
# cannot be an assay difference. The cost is real and should be stated when the results are:
# HNRNPC-Hela-iCLIP was the only assay-matched column (THRAP3 is itself iCLIP), and HNRNPC
# agreeing with ITSELF across both cell line and assay was the strongest result of the earlier
# runs. HNRNPC keeps its HepG2 and K562 columns, but that replication is now within-assay.
# TRA2B-MDA231-iCLIP is dropped too. Point -s back at RBPeekSamplesheet.tsv to restore both.
#
# --no-clustering. k-means on binarised presence answers "what kinds of locus are there",
# which is a Centrosome question; here the question is which panel columns are enriched at
# THRAP3 sites. Dropping it removes binf_heatmap_clusters.tsv and the six per-cluster
# metaprofiles, orders the heatmap's locus axis by total support instead, and leaves the tSNE
# a single colour. protein_ranking.tsv is the replacement output: every panel sample ranked by
# mean peak support, carrying frac_centred alongside, so depth and central enrichment can be
# compared directly rather than one being chosen in advance.
#
# No -i and no --protein-nt-heatmap. The metaprofile now takes its top 10 from the SAME
# enrichment ranking that picks the heatmap's 20 columns, so it already shows the heatmap's
# proteins as profiles - which was the only job protein_nt_metaprofile_heatmap had. --window
# stays at 100 throughout; nothing about the statistics changes.
# --metaprofile-top-proteins 20 matches --cluster-top-proteins 20, so the metaprofile and the
# heatmap now show the IDENTICAL set rather than the metaprofile showing the top half of it.
# At 10 they diverged in a way that read as a bug: enrichment rank and depth are only weakly
# related here (Spearman +0.35 over the 224 eligible columns), so the metaprofile's ranks 1-10
# are mostly shallow columns while the heatmap's brightest rows - BCLAF1, BUD13, SAFB - sit at
# enrichment ranks 16, 17 and 12 and were drawn on the heatmap but absent from the curves.
python3 scripts/intersect_inference_bed.py \
  -x ../CLIP \
  -b "$BINF" \
  -s THRAP3/RBPeekSamplesheet_eCLIP.tsv \
  --skip-merge \
  --table \
  -o "$OUT" \
  --window 100 \
  --panel-anchor midpoint \
  --protein-select enrichment \
  --enrichment-window 5 \
  --centrality-min-total 100 \
  --protein-min-loci-frac 0.0248 \
  --heatmap-min-support-percentile 10 \
  --heatmap-scale percentile \
  --heatmap-scale-percentile 99 \
  --no-clustering \
  --cluster-top-proteins 20 \
  --metaprofile-top-proteins 20 \
  --tsne
