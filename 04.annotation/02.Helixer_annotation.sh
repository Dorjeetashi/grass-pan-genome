#!/usr/bin/env bash
set -euo pipefail

# Path to the Singularity container (modify according to your actual path)
SIF="~/singularity/helixer_v0.3.5.sif"

# Paths to the three species' common path FASTA files
# Replace these with your actual file paths
FASTA_SOR="SOR.common.fa"          # Sorghum common path genome
FASTA_RICE="RICE.common.fa"        # Rice common path genome
FASTA_MAIZE="MAIZE.common.fa"      # Maize common path genome

# Output directory for Helixer annotations (recommended)
OUT_DIR="./helixer_annotations"
mkdir -p "${OUT_DIR}"

# Common parameters
LINEAGE="land_plant"
MODEL_DIR="~/.local/share/Helixer/models/"
TMP_BASE_DIR="./helixer_tmp"
mkdir -p "${TMP_BASE_DIR}"

# ----------------------------------------------------------------------
# 1. Annotate Sorghum common path genome
# ----------------------------------------------------------------------
echo "Annotating sorghum common path genome..."
singularity exec "${SIF}" Helixer.py \
  --lineage "${LINEAGE}" \
  --downloaded-model-path "${MODEL_DIR}" \
  --fasta-path "${FASTA_SOR}" \
  --gff-output-path "${OUT_DIR}/sorghum.common.helixer.gff3" \
  --temporary-dir "${TMP_BASE_DIR}/sorghum" \
  --species sorghum_common

# ----------------------------------------------------------------------
# 2. Annotate Rice common path genome
# ----------------------------------------------------------------------
echo "Annotating rice common path genome..."
singularity exec "${SIF}" Helixer.py \
  --lineage "${LINEAGE}" \
  --downloaded-model-path "${MODEL_DIR}" \
  --fasta-path "${FASTA_RICE}" \
  --gff-output-path "${OUT_DIR}/rice.common.helixer.gff3" \
  --temporary-dir "${TMP_BASE_DIR}/rice" \
  --species rice_common

# ----------------------------------------------------------------------
# 3. Annotate Maize common path genome
# ----------------------------------------------------------------------
echo "Annotating maize common path genome..."
singularity exec "${SIF}" Helixer.py \
  --lineage "${LINEAGE}" \
  --downloaded-model-path "${MODEL_DIR}" \
  --fasta-path "${FASTA_MAIZE}" \
  --gff-output-path "${OUT_DIR}/maize.common.helixer.gff3" \
  --temporary-dir "${TMP_BASE_DIR}/maize" \
  --species maize_common

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
echo ""
echo "All Helixer annotations completed."
echo "Output files:"
ls -lh "${OUT_DIR}"/*.gff3
echo ""
echo "Temporary directories (can be deleted after checking results):"
ls -ld "${TMP_BASE_DIR}"/*
