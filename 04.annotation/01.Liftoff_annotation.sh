#!/usr/bin/env bash
set -euo pipefail

# Container path
SIF="~/singularity/GENETools202309.sif"

# Array of species (short name, fasta, annotation)
declare -A species=(
    ["sorghum"]="SOR.common.fa sorghum_genes.gff3"
    ["rice"]="RICE.common.fa rice_genes.gff3"
    ["maize"]="MAIZE.common.fa maize_genes.gff3"
)

for sp in "${!species[@]}"; do
    read fasta annot <<< "${species[$sp]}"
    
    prefix="${sp}.common"
    
    echo "Running Liftoff for ${sp^} common path genome..."
    
    singularity exec "${SIF}" liftoff \
        -polish \
        -u "${prefix}.unmapped.txt" \
        -g "${prefix}.liftoff.gff" \
        -o "${prefix}.liftoff.gff" \
        "${fasta}" \
        "${annot}"
done

echo ""
echo "Liftoff finished for all species."
ls -lh *.common.liftoff.gff *.unmapped.txt
