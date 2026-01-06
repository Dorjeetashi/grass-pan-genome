# Soft-mask genomes
repeatMasker -e ncbi -species "sorghum bicolor" -gff -xsmall -html -norna -source -no_is -pa 40 -s Sb_common.fasta
repeatMasker -e ncbi -species rice -gff -xsmall -html -norna -source -no_is -pa 40 -s Os_common.fasta
repeatMasker -e ncbi -species "Zea mays" -gff -xsmall -html -norna -source -no_is -pa 40 -s Zm_common.fasta
# Generate BAM files for Sorghum RNA-seq data
hisat2-build Sb.common.sm.fa Sb.common
for r1 in /public/home/01.data/03.pan_genome/annotation/new_common/Sb_rna/Sb_RNAseq/*.1.clean.fq.gz; do
    # Extract sample ID from filename
    base=$(basename "$r1" .1.clean.fq.gz)
    r2=/public/home/01.data/03.pan_genome/annotation/new_common/Sb_rna/Sb_RNAseq/${base}.2.clean.fq.gz

    # Check if the corresponding R2 file exists
    if [[ -f "$r2" ]]; then
        echo "Processing $base..."

        # Run hisat2 and samtools sort
        hisat2 --dta --new-summary -p 40 -x Sb.common -1 "$r1" -2 "$r2" 2> "${base}.hisat2.log" | samtools sort -@ 10 > "${base}.bam"

        # Index BAM file
        samtools index "${base}.bam"
    else
        echo "Warning: $r2 not found, skipping $base"
    fi
done
## Sorghum common annotation with BRAKER

# First, create a file containing all Sorghum BAM paths (once)
ls /public/home/01.data/03.pan_genome/annotation/new_common/Sb_rna/Sb_common_bam/common_bam/*.bam \
  > Sb_common_bams.list

# Join them into a comma-separated list for BRAKER
sb_bam_list=$(paste -sd, Sb_common_bams.list)

singularity exec /public/home/singularity_package/braker3.sif braker.pl \
  --threads=40 \
  --species=BTx623 \
  --AUGUSTUS_CONFIG_PATH=/public/home/01.data/03.pan_genome/annotation/new_common/augustus_config \
  --genome=Sb.common.sm.fa \
  --prot_seq=/public/home/01.data/03.pan_genome/annotation/new_common/Os_rna/Sb_cluster.fasta \
  --useexisting \
  --bam="${sb_bam_list}"

# Generate BAM files for Rice RNA-seq data
hisat2-build Sb.common.sm.fa Os.common
for r1 in /public/home/01.data/03.pan_genome/annotation/new_common/Os_rna/Os_RNAseq/CRR*_f1.fq.gz; do
    # Extract sample ID from filename
    base=$(basename "$r1" _f1.fq.gz)
    r2=/public/home/01.data/03.pan_genome/annotation/new_common/Os_rna/Os_RNAseq/${base}_r2.fq.gz

    # Check if the corresponding R2 file exists
    if [[ -f "$r2" ]]; then
        echo "Processing $base..."

        # Run hisat2 and samtools sort
        hisat2 --dta --new-summary -p 40 -x Os.common -1 "$r1" -2 "$r2" 2> "${base}.common.hisat2.log" | samtools sort -@ 20 > "${base}.common.bam"

        # Index BAM file
        samtools index "${base}.common.bam"
    else
        echo "Warning: $r2 not found, skipping $base"
    fi
done
## Rice common annotation with BRAKER

ls /public/home/01.data/03.pan_genome/annotation/new_common/Os_rna/Os_common_bam/*.bam \
  > Os_common_bams.list

os_bam_list=$(paste -sd, Os_common_bams.list)

singularity exec /public/home/singularity_package/braker3.sif braker.pl \
  --threads=40 \
  --species=rice \
  --AUGUSTUS_CONFIG_PATH=/public/home/01.data/03.pan_genome/annotation/new_common/config \
  --genome=Os.common.sm.fa \
  --prot_seq=/public/home/01.data/03.pan_genome/annotation/new_common/Os_rna/Os_cluster.fasta \
  --useexisting \
  --softmasking \
  --bam="${os_bam_list}"

# Generate BAM files for Maize RNA-seq data
hisat2-build Zm.common.sm.fa Zm.common
for r1 in /public/home/01.data/03.pan_genome/annotation/new_common/Zm_rna/fastq_output/RNA_seq/SRR*_1.fastq.gz; do
    # Extract sample ID from filename
    base=$(basename "$r1" _1.fastq.gz)
    r2=/public/home/01.data/03.pan_genome/annotation/new_common/Zm_rna/fastq_output/RNA_seq/${base}_2.fastq.gz

    # Check if the corresponding R2 file exists
    if [[ -f "$r2" ]]; then
        echo "Processing $base..."

        # Run hisat2 and samtools sort
        hisat2 --dta --new-summary -p 40 -x Zm.common -1 "$r1" -2 "$r2" 2> "${base}.hisat2.log" | samtools sort -@ 20 > "${base}.bam"

        # Index BAM file
        samtools index "${base}.bam"
    else
        echo "Warning: $r2 not found, skipping $base"
    fi
done
## Maize common annotation with BRAKER (cleaned version)

ls /public/home/01.data/03.pan_genome/annotation/new_common/Zm_rna/fastq_output/Zm_common/*.bam \
  > Zm_common_bams.list

zm_bam_list=$(paste -sd, Zm_common_bams.list)

singularity exec /public/home/singularity_package/braker3.sif braker.pl \
  --threads=40 \
  --species=maize \
  --AUGUSTUS_CONFIG_PATH=/public/home/01.data/03.pan_genome/annotation/new_common/config_new \
  --skipAllTraining \
  --genome=Zm.common.sm.fa \
  --useexisting \
  --softmasking \
  --bam="${zm_bam_list}"
