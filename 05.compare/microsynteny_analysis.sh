#!/usr/bin/env bash
set -euo pipefail

########################################
# Step 1: Convert GFF3 annotations to BED
########################################

# Format: "GFF3_filename  BED_prefix"
species_list=(
  "Sorghum_BTx623.longest_isoform.gff3 Sb_BTx623"
  "Sorghum_IS1661.longest_isoform.gff3 Sb_IS1661"
  "Sorghum_Ji2731.longest_isoform.gff3 Sb_Ji2731"
  "Sorghum_Tx430.longest_isoform.gff3 Sb_Tx430"
  "Sorghum_PI525695.longest_isoform.gff3 Sb_PI525695"
  "Zm_B73.longest_isoform.gff3 Zm_B73"
  "Os_MSU.longest_isoform.gff3 Os_MSU"
)

echo ">>> Converting GFF3 to BED ..."
for item in "${species_list[@]}"; do
  gff3_file=$(echo "$item" | awk '{print $1}')
  prefix=$(echo "$item" | awk '{print $2}')
  echo "  - ${gff3_file} -> ${prefix}.bed"
  python3 -m jcvi.formats.gff bed \
    --type=mRNA \
    --key=ID \
    "${gff3_file}" \
    -o "${prefix}.bed"
done

########################################
# Step 2: Merge selected BED files
########################################

echo ">>> Merging BED files into all.bed ..."
python3 -m jcvi.formats.bed merge \
  Sb_BTx623.bed \
  Sb_IS1661.bed \
  Sb_Ji2731.bed \
  Sb_Tx430.bed \
  Sb_PI525695.bed \
  Zm_B73.bed \
  -o all.bed

########################################
# Step 3: Ortholog detection (Sb_BTx623 as reference)
########################################

ref="Sb_BTx623"
targets=(
  "Sb_IS1661"
  "Sb_Ji2731"
  "Sb_Tx430"
  "Sb_PI525695"
  "Zm_B73"
  "Os_MSU"
)

echo ">>> Running jcvi.compara.catalog ortholog (reference = ${ref}) ..."
for tgt in "${targets[@]}"; do
  echo "  - ${ref} vs ${tgt}"
  python3 -m jcvi.compara.catalog ortholog \
    "${ref}" "${tgt}" \
    --no_strip_names
done

########################################
# Step 4: Identify syntenic blocks by MCscan
########################################

echo ">>> Running jcvi.compara.synteny mcscan ..."
for tgt in "${targets[@]}"; do
  anchors_file="${ref}.${tgt}.anchors"
  block_file="${ref}.${tgt}.block"
  echo "  - ${anchors_file} -> ${block_file}"
  python3 -m jcvi.compara.synteny mcscan \
    "${ref}.bed" "${anchors_file}" \
    --iter=1 \
    -o "${block_file}"
done

########################################
# Step 5: Merge block files and extract SbOMT3 region
########################################

echo ">>> Joining block files into all.blocks ..."
python3 -m jcvi.formats.base join \
  Sb_BTx623.Sb_IS1661.block \
  Sb_BTx623.Sb_Ji2731.block \
  Sb_BTx623.Sb_Tx430.block \
  Sb_BTx623.Sb_PI525695.block \
  Sb_BTx623.Zm_B73.block \
  Sb_BTx623.Os_MSU.block \
  --noheader | cut -f1,2,4,6,8,10,12 > all.blocks

echo ">>> Extracting syntenic region containing SbOMT3 ..."
grep -10 "Sobic.006G007900.1.v5.1" all.blocks > SbOMT3_block

########################################
# Step 6: Prepare layout file and plot
########################################

echo ">>> Writing layout file SbOMT3.layout ..."
cat > SbOMT3.layout << 'EOF'
#x, y, rotation, ha, va, color, ratio, label
0.5, 0.9, 0, left, center, pink, 1.5, Sb_BTx623
0.5, 0.8, 0, left, center, pink, 1.5, Sb_IS1661
0.5, 0.7, 0, left, center, pink, 1.8, Sb_Ji2731
0.5, 0.6, 0, left, center, pink, 1.3, Sb_Tx430
0.5, 0.5, 0, left, center, pink, 2.5, Sb_PI525695
0.5, 0.4, 0, left, center, pink, .6, Zm_B73
0.5, 0.3, 0, left, center, pink, 4.2, Os_MSU
#edges
e, 0, 1
e, 1, 2
e, 2, 3
e, 3, 4
e, 4, 5
e, 5, 6
EOF

echo ">>> Plotting synteny figure ..."
python3 -m jcvi.graphics.synteny \
  SbOMT3_block \
  all.bed \
  SbOMT3.layout \
  --genelabelsize=6 \
  --genelabels=OMT3 \
  --shadestyle=line \
  --glyphcolor=orthogroup

echo ">>> Done."
