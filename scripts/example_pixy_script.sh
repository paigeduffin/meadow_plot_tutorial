#!/bin/bash
#SBATCH --job-name=pixy_windowed_pi_fst
#SBATCH --partition=main
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=pixy_windowed_pi_fst.%j.out
#SBATCH --error=pixy_windowed_pi_fst.%j.err

# This template was written by Paige Duffin (2026) to accompany the meadow 
# plot tutorial available at: https://paigeduffin.github.io/meadow_plot_tutorial
# 
# this template uses pixy, a valuable tool external to the tutorial: 
# Korunes, K. L. & Samuk, K. pixy: Unbiased estimation of nucleotide diversity 
# and divergence in the presence of missing data. Mol. Ecol. Resour. 21, 
# 1359–1368 (2021). https://doi.org/10.1111/1755-0998.13326
#
# This template calculates windowed nucleotide diversity (pi) and Fst with pixy.
# Lines that must be reviewed or edited by the user are marked with "EDITME".
# Slurm log files are written to the directory from which this script is submitted.

# See pixy documentation for important additional info needed to run the program 
# at: https://pixy.readthedocs.io/en/latest/index.html

###############################################################################
# Preparing the sites and windows files
#
# SITES FILE (optional)
# A sites file restricts the analysis to specified genomic positions. It must
# be headerless and contain two tab-separated columns: chromosome and position.
# To create one from a filtered VCF containing the sites you want retained:
# 
# bcftools query -f '%CHROM\t%POS\n' filtered_sites.vcf.gz > sites.txt
#
# As an example, Duffin, Ruggeri et al. (in prep.) generated a sites file from
# a VCF previously filtered to exclude SNPs that violated Hardy–Weinberg
# equilibrium (HWE). If every site in VCF_FILE should be analyzed, a sites file
# is unnecessary. Remove SITES_FILE and the --sites_file argument below.
#
# WINDOWS (BED) FILE
# First index the reference genome and create a chromosome-length file:
#
# samtools faidx reference.fasta
# cut -f1,2 reference.fasta.fai > genome.sizes
#
# Then create a BED file specifying the desired windows. For example, Duffin,
# Using sites.txt as the list of callable positions, create a BED file specifying
# the desired windows. For example, Duffin, Ruggeri et al. (in prep.) created
# nonoverlapping 10-kb windows and retained only windows containing at least
# 5 kb of callable sites (50%):
#
# awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' sites.txt > callable_sites.bed
# bedtools makewindows -g genome.sizes -w 10000 |
#   bedtools intersect -a - -b callable_sites.bed -c |
#   awk 'BEGIN{OFS="\t"} $4 >= 5000 {print $1,$2,$3}' > windows_10kb_min5kb_callable.bed
#
# Add -s #### after -w 10000 to specify a step size and create overlapping windows.
#
# Move the resulting files into the input directory specified below.
###############################################################################

set -euo pipefail

###############################################################################
# User settings

PROJECT_DIR="/path/to/meadow_plot_project"    # EDITME: main project directory
OUTPUT_DIR="${PROJECT_DIR}/pixy"
POPULATIONS_FILE="${PROJECT_DIR}/input/populations.txt"  # EDITME if the filename differs

VCF_FILE="/path/to/all_sites.vcf.gz"          # EDITME: bgzipped and indexed all-sites VCF
SITES_FILE="/path/to/sites.txt"               # EDITME: tab-separated CHROM and POS file

# This BED file defines the analyzed windows. 
# The meadow-plot tutorial uses 15-kb windows with a 5-kb step.
WINDOWS_FILE="${PROJECT_DIR}/input/windows_15kb_5kb_step.bed"  # EDITME if the filename differs

CONDA_MODULE="conda/25.11.0"                  # EDITME: module available on your cluster
PIXY_ENV="/path/to/conda/envs/pixy"           # EDITME: pixy conda environment
OUTPUT_PREFIX="15kb_5kb_step"
N_CORES="${SLURM_CPUS_PER_TASK:-4}"

###############################################################################
# Prepare output directory

mkdir -p "$OUTPUT_DIR"

###############################################################################
# Load and activate pixy

module load "$CONDA_MODULE"

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$PIXY_ENV"

###############################################################################
# Calculate pi and Fst

pixy --stats pi fst \
  --vcf "$VCF_FILE" \
  --populations "$POPULATIONS_FILE" \
  --include_multiallelic_snps \
  --sites_file "$SITES_FILE" \
  --bed_file "$WINDOWS_FILE" \
  --n_cores "$N_CORES" \
  --output_prefix "$OUTPUT_PREFIX" \
  --output_folder "$OUTPUT_DIR"
