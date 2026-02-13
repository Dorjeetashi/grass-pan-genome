# Building "Common Path" Consensus Genomes from Pangenome VCFs

This document describes the steps to:

1. Prepare reference genomes and pangenome VCFs for sorghum (Sb), rice (Os), and maize (Zm).
2. Normalize multi-allelic VCFs.
3. Generate a "common path" single-sample VCF by selecting major alleles (present in >50% of samples, SVLEN ≥ 500 bp by default).
4. Create consensus FASTA sequences representing the "common path" genome.
5. Identify structural variants (SVs) that overlap genes.

## 1. Software Installation

```bash
conda create -n common -c bioconda pysam bcftools samtools tqdm bedtools
conda activate common
#2. Input Files 

#3.1 Index reference genomes
# Sorghum
REF_SOR=/path/to/sorghum.ref.fa
VCF_SOR=/path/to/sorghum.pangenome.vcf.gz

# Rice
REF_RICE=/path/to/rice.ref.fa
VCF_RICE=/path/to/rice.pangenome.vcf.gz

# Maize
REF_MAIZE=/path/to/maize.ref.fa
VCF_MAIZE=/path/to/maize.pangenome.vcf.gz
samtools faidx "${REF_SOR}"
samtools faidx "${REF_RICE}"
samtools faidx "${REF_MAIZE}"

#3.2 Normalize VCFs (left-align, split multi-allelic sites, add statistics)
for SPECIES in SOR RICE MAIZE; do
  REF=$(eval echo \$REF_${SPECIES})
  VCF=$(eval echo \$VCF_${SPECIES})

  bcftools norm -f "${REF}" -m -any -Oz -o "${SPECIES}.norm.vcf.gz" "${VCF}"
  bcftools index -t "${SPECIES}.norm.vcf.gz"

  # Sort (optional but recommended)
  zcat "${SPECIES}.norm.vcf.gz" | bcftools sort -Oz -o "${SPECIES}.norm.sorted.vcf.gz"
  bcftools index -t "${SPECIES}.norm.sorted.vcf.gz"
done

#3.3Create make_common_path.py script
chmod +x make_common_path.py
3.4 Run the script for each species
./make_common_path.py SOR.norm.sorted.vcf.gz "${REF_SOR}" SOR.common.vcf.gz SOR.common.provenance.tsv 500 0.5
./make_common_path.py RICE.norm.sorted.vcf.gz "${REF_RICE}" RICE.common.vcf.gz RICE.common.provenance.tsv 500 0.5
./make_common_path.py MAIZE.norm.sorted.vcf.gz "${REF_MAIZE}" MAIZE.common.vcf.gz MAIZE.common.provenance.tsv 500 0.5

#3.5 Generate consensus "common path" FASTA files
for SPECIES in SOR RICE MAIZE; do
  case ${SPECIES} in
    SOR) REF="${REF_SOR}" ;;
    RICE) REF="${REF_RICE}" ;;
    MAIZE) REF="${REF_MAIZE}" ;;
  esac

  VCF="${SPECIES}.common.vcf.gz"

  # Index the new VCF
  bcftools index -t "${VCF}"

  # Generate consensus sequence (applies GT=1/1 for selected ALTs)
  bcftools consensus \
    -f "${REF}" \
    -s COMMON \
    -H 1 \
    "${VCF}" \
    > "${SPECIES}.common.fa"

  # Optional: index the FASTA
  samtools faidx "${SPECIES}.common.fa"
done

#4. Quick inspection: Which genes are affected by selected SVs?
# Convert provenance TSV → BED (chrom, start, end, svtype, svlen)
awk 'NR>1 {print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' SOR.common.provenance.tsv > SOR.common_sv.bed
awk 'NR>1 {print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' RICE.common.provenance.tsv > RICE.common_sv.bed
awk 'NR>1 {print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' MAIZE.common.provenance.tsv > MAIZE.common_sv.bed

# Example: intersect with gene annotations (adjust GFF path)
grep -w "gene" MSU_all_genes.gff3 > MSU_genes.gff3

bedtools intersect -a SOR.common_sv.bed -b sorghum_genes.gff3 -wa -wb > SOR_sv_gene_overlap.tsv
bedtools intersect -a RICE.common_sv.bed -b rice_genes.gff3 -wa -wb > RICE_sv_gene_overlap.tsv
bedtools intersect -a MAIZE.common_sv.bed -b maize_genes.gff3 -wa -wb > MAIZE_sv_gene_overlap.tsv


