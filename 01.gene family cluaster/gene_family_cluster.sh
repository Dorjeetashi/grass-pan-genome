### Data Preparation
ln -s /public/home/zhaxi/grass_pangenome/01.ortho_ref/grass_pep/OrthoFinder/Results_May28/Orthogroups/Orthogroups.GeneCount.tsv ./ref.GeneCount.tsv

## Single reference genome Venn diagram
sed '1d' ref.GeneCount.tsv | awk '{if($2 > 0 ) print $1}' > Os.id
sed '1d' ref.GeneCount.tsv | awk '{if($3 > 0 ) print $1}' > Sb.id
sed '1d' ref.GeneCount.tsv | awk '{if($4 > 0 ) print $1}' > Zm.id

## Plotting:
http://bioinformatics.psb.ugent.be/webtools/Venn/

## Pan-genome Venn diagram
ln -s /public/home/zhaxi/grass_pangenome/ortho_pan/grass_pep_pan/OrthoFinder/Results_June10/Orthogroups/Orthogroups.GeneCount.tsv ./pan.GeneCount.tsv
# Extract OG IDs for species starting with Sb_

awk 'NR==1 {for(i=2;i<=NF;i++) if($i~/^Sb_/) sb[i]=1} 
     NR>1 {for(i in sb) if($i>=1) {print $1; next}}' pan.GeneCount.tsv > Sb_genefamily.ID.txt

# Extract OG IDs for species starting with Os_
awk 'NR==1 {for(i=2;i<=NF;i++) if($i~/^Os_/) os[i]=1} 
     NR>1 {for(i in os) if($i>=1) {print $1; next}}' pan.GeneCount.tsv > Os_genefamily.ID.txt

# Extract OG IDs for species starting with Zm_
awk 'NR==1 {for(i=2;i<=NF;i++) if($i~/^Zm_/) zm[i]=1} 
     NR>1 {for(i in zm) if($i>=1) {print $1; next}}' pan.GeneCount.tsv > Zm_genefamily.ID.txt

## Plotting:
http://bioinformatics.psb.ugent.be/webtools/Venn/

## Extract common species genes
comm -12 Sb_genefamily.ID.txt Os_genefamily.ID.txt | comm -12 - Zm_genefamily.ID.txt > common_all.txt

## Extract common genes between Sorghum and Maize
comm -12 Sb_genefamily.ID.txt Zm_genefamily.ID.txt | comm -23 - common_all.txt > Sb_Zm_common.txt

## Extract common genes between Sorghum and Rice
comm -12 Sb_genefamily.ID.txt Os_genefamily.ID.txt | comm -23 - common_all.txt > Sb_Os_common.txt

## Extract common genes between Maize and Rice
comm -12 Zm_genefamily.ID.txt Os_genefamily.ID.txt | comm -23 - common_all.txt > Zm_Os_common.txt

## Extract species-specific genes
# Sorghum-specific IDs (present in Sorghum but not in Maize or Rice)
grep -vx -f <(cat Zm_genefamily.ID.txt Os_genefamily.ID.txt | sort -u) Sb_genefamily.ID.txt > Sb_genefamily.unique.txt

# Maize-specific IDs (present in Maize but not in Sorghum or Rice)
grep -vx -f <(cat Sb_genefamily.ID.txt Os_genefamily.ID.txt | sort -u) Zm_genefamily.ID.txt > Zm_genefamily.unique.txt

# Rice-specific IDs (present in Rice but not in Maize or Sorghum)
grep -vx -f <(cat Sb_genefamily.ID.txt Zm_genefamily.ID.txt | sort -u) Os_genefamily.ID.txt > Os_genefamily.unique.txt

## Gene family classification in a single species
## In Sorghum, present in at least 16 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Sb/) col[i]=1; totalSb=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= totalSb - 0) print $1}' Orthogroups.GeneCount.tsv > Sb.core_gene_families.txt

## In Sorghum, present in 2-15 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Sb/) col[i]=1; totalSb=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= 2 && count <= 15) print $1}' Orthogroups.GeneCount.tsv > Sb.shell_gene_families.txt

## In Sorghum, present in only 1 variety
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Sb/) col[i]=1; totalSb=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count == 1) print $1}' Orthogroups.GeneCount.tsv > Sb.cloud_gene_families.txt

## In Maize, present in 26 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Zm/) col[i]=1; totalZm=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= totalZm - 0) print $1}' Orthogroups.GeneCount.tsv > Zm.core_gene_families.txt

## In Maize, present in 2-25 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Zm/) col[i]=1; totalZm=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= 2 && count <= 25) print $1}' Orthogroups.GeneCount.tsv > Zm.shell_gene_families.txt

## In Maize, present in only 1 variety
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Zm/) col[i]=1; totalZm=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count == 1) print $1}' Orthogroups.GeneCount.tsv > Zm.cloud_gene_families.txt

## In Rice, present in 33 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Os/) col[i]=1; totalOs=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= totalOs - 0) print $1}' Orthogroups.GeneCount.tsv > Os.core_gene_families.txt

## In Rice, present in 2-32 varieties
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Os/) col[i]=1; totalOs=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count >= 2 && count <= 32) print $1}' Orthogroups.GeneCount.tsv > Os.shell_gene_families.txt

## In Rice, present in only 1 variety
awk 'NR==1{for(i=2;i<=NF;i++) if($i ~ /^Os/) col[i]=1; totalOs=length(col)} NR>1{count=0; for(i in col) if($i > 0) count++; if(count == 1) print $1}' Orthogroups.GeneCount.tsv > Os.cloud_gene_families.txt

## Extract core common gene families
comm -12 Sb.core_gene_families.txt Zm.core_gene_families.txt | comm -12 - Os.core_gene_families.txt > core.common_all.txt

## Extract core common gene families between Sorghum and Maize
comm -12 Sb.core_gene_families.txt Zm.core_gene_families.txt | comm -23 - core.common_all.txt > Sb_Zm_core.common.txt

## Extract core common gene families between Sorghum and Rice
comm -12 Sb.core_gene_families.txt Os.core_gene_families.txt | comm -23 - core.common_all.txt > Sb_Os_core.common.txt

## Extract core common gene families between Maize and Rice
comm -12 Zm.core_gene_families.txt Os.core_gene_families.txt | comm -23 - core.common_all.txt > Zm_Os_core.common.txt

## Extract species-specific core gene families
# Sorghum-specific OG IDs (present in Sorghum but not in Maize or Rice)
grep -vx -f <(cat Zm.core_gene_families.txt Os.core_gene_families.txt | sort -u) Sb.core_gene_families.txt > Sb_core.unique.txt

# Maize-specific OG IDs
grep -vx -f <(cat Sb.core_gene_families.txt Os.core_gene_families.txt | sort -u) Zm.core_gene_families.txt > Zm_core.unique.txt

# Rice-specific OG IDs
grep -vx -f <(cat Sb.core_gene_families.txt Zm.core_gene_families.txt | sort -u) Os.core_gene_families.txt > Os_core.unique.txt

## Extract shell common gene families 
comm -12 Sb.shell_gene_families.txt Zm.shell_gene_families.txt | comm -12 - Os.shell_gene_families.txt > shell.common_all.txt

## Extract shell common gene families between Sorghum and Maize
comm -12 Sb.shell_gene_families.txt Zm.shell_gene_families.txt | comm -23 - shell.common_all.txt > Sb_Zm_shell.common.txt

## Extract shell common gene families between Sorghum and Rice
comm -12 Sb.shell_gene_families.txt Os.shell_gene_families.txt | comm -23 - shell.common_all.txt > Sb_Os_shell.common.txt

## Extract shell common gene families between Maize and Rice
comm -12 Zm.shell_gene_families.txt Os.shell_gene_families.txt | comm -23 - shell.common_all.txt > Zm_Os_shell.common.txt

## Extract species-specific shell gene families
# Sorghum-specific OG IDs (present in Sorghum but not in Maize or Rice)
grep -vx -f <(cat Zm.shell_gene_families.txt Os.shell_gene_families.txt | sort -u) Sb.shell_gene_families.txt > Sb_shell.unique.txt

# Maize-specific OG IDs
grep -vx -f <(cat Sb.shell_gene_families.txt Os.shell_gene_families.txt | sort -u) Zm.shell_gene_families.txt > Zm_shell.unique.txt

# Rice-specific OG IDs
grep -vx -f <(cat Sb.shell_gene_families.txt Zm.shell_gene_families.txt | sort -u) Os.shell_gene_families.txt > Os_shell.unique.txt

## Extract cloud common gene families
comm -12 Sb.cloud_gene_families.txt Zm.cloud_gene_families.txt | comm -12 - Os.cloud_gene_families.txt > cloud.common_all.txt

## Extract cloud common gene families between Sorghum and Maize
comm -12 Sb.cloud_gene_families.txt Zm.cloud_gene_families.txt | comm -23 - cloud.common_all.txt > Sb_Zm_cloud.common.txt

## Extract cloud common gene families between Sorghum and Rice
comm -12 Sb.cloud_gene_families.txt Os.cloud_gene_families.txt | comm -23 - cloud.common_all.txt > Sb_Os_cloud.common.txt

## Extract cloud common gene families between Maize and Rice
comm -12 Zm.cloud_gene_families.txt Os.cloud_gene_families.txt | comm -23 - cloud.common_all.txt > Zm_Os_cloud.common.txt

## Extract species-specific cloud gene families
# Sorghum-specific OG IDs (present in Sorghum but not in Maize or Rice)
grep -vx -f <(cat Zm.cloud_gene_families.txt Os.cloud_gene_families.txt | sort -u) Sb.cloud_gene_families.txt > Sb_cloud.unique.txt

# Maize-specific OG IDs
grep -vx -f <(cat Sb.cloud_gene_families.txt Os.cloud_gene_families.txt | sort -u) Zm.cloud_gene_families.txt > Zm_cloud.unique.txt

# Rice-specific OG IDs
grep -vx -f <(cat Sb.cloud_gene_families.txt Zm.cloud_gene_families.txt | sort -u) Os.cloud_gene_families.txt > Os_cloud.unique.txt

